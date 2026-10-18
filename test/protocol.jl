include("protocol_fixtures.jl")
using .ProtocolFixtures: DerivativeProbe, HistoryProbe, PreparationCounts, CountedReal, consume!

# Independent callback observations, including fault injection at a chosen call.
# These allocate deliberately; they are behavioral tests, not performance data.
function observed_callbacks(T; model=:linear)
    calls = Dict(:values=>0,:derivative=>0,:joint=>0)
    events,points = Symbol[],Vector{T}[]
    fault = Ref((operation=:none,call=0,action=:ok))
    values(F,y) = fill!(F,model===:linear ? -sum(y) : sum(v -> (v-one(T))^2,y))
    derivatives(D,y,d) = fill!(D,model===:linear ? -sum(d) :
        T(2)*sum((y[i]-one(T))*d[i] for i in eachindex(y)))
    function observe(operation,y,F,D)
        calls[operation]+=1
        push!(events,operation); push!(points,collect(y))
        f=fault[]
        if operation===f.operation && calls[operation]==f.call
            f.action===:error && return SubString("!reported failure!",2,17)
            f.action===:throw && error("callback exception")
            f.action===:invalid && return 42
            f.action===:nonfinite_values && (F[end]=T(NaN))
            f.action===:nonfinite_derivatives && (D[end]=T(Inf))
        end
        return nothing
    end
    value! = (F,y) -> (values(F,y); observe(:values,y,F,nothing))
    derivative! = (D,y,d) -> (derivatives(D,y,d); observe(:derivative,y,nothing,D))
    joint! = (F,D,y,d) -> (values(F,y); derivatives(D,y,d); observe(:joint,y,F,D))
    return (;value!,derivative!,joint!,calls,events,points,fault)
end

counts(r) = (r.nfev,r.ndev)
observed_counts(cb) = (cb.calls[:values]+cb.calls[:joint], cb.calls[:derivative]+cb.calls[:joint])
scalars(r) = (r.success,r.alpha,r.reason,r.message,counts(r),r.iterations)
copy_accept!(x,Fx,r) = (copyto!(x,r.x_trial);copyto!(Fx,r.F_trial);nothing)

@testset "One consumer, three methods ($T)" for T in (Float32,Float64)
    for spec in (BackTracking(),DerivativeProbe(),HistoryProbe())
        cb=observed_callbacks(T)
        x,Fx = zeros(T,2),zeros(T,3)
        search=configure_linesearch(spec,cb.value!,x,Fx;
            derivative! = cb.derivative!,value_derivative! = cb.joint!)
        @test isempty(cb.events)
        steps=[(fill(T(a),2),spec isa BackTracking ?
                (slopes=fill(-T(2)*T(a),3),factor=-T(2)*T(a)) : (;))
               for a in (1,0.5,0.25)]
        results=consume!(search,x,Fx,steps)
        @test length(results)==3 && all(r -> r.success,results)
        @test x==fill(T(1.75),2) && Fx==fill(T(-3.5),3)
        @test all(r -> r.x_trial===search.x_trial && r.F_trial===search.F_trial,results)
        @test results[1].x_trial==x && results[1].F_trial==Fx # borrowed, not snapshots
        if spec isa DerivativeProbe
            @test all(r -> counts(r)==(1,2),results)
            @test cb.events==repeat([:derivative,:joint],3)
        else
            @test all(r -> counts(r)==(1,0),results)
            @test cb.calls[:derivative]==cb.calls[:joint]==0
        end
        @test (sum(r.nfev for r in results),sum(r.ndev for r in results))==observed_counts(cb)
        # A failure after accepted/copied successes leaves earlier scalars intact.
        first_scalars=scalars(results[1]); accepted=(copy(x),copy(Fx))
        cb.fault[]=(operation=spec isa DerivativeProbe ? :joint : :values,
            call=spec isa DerivativeProbe ? cb.calls[:joint]+1 : cb.calls[:values]+1,
            action=:error)
        failed=only(consume!(search,x,Fx,[(ones(T,2),steps[1][2])];maxstep=0.5))
        @test !failed.success && failed.reason===:evaluation_failure
        @test (x,Fx)==accepted && scalars(results[1])==first_scalars
        @test failed.x_trial===results[1].x_trial && failed.F_trial===results[1].F_trial
        @test failed.alpha==T(0.5) && results[1].alpha==one(T)
    end
end

@testset "Derivative capabilities and exact work ($T)" for T in (Float32,Float64)
    x,Fx,d=zeros(T,1),zeros(T,2),ones(T,1)
    cb=observed_callbacks(T)
    @test_throws ArgumentError configure_linesearch(DerivativeProbe(),cb.value!,x,Fx)
    @test_throws ArgumentError configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;derivative! = identity)
    @test_throws ArgumentError configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;value_derivative! = identity)
    @test_throws ArgumentError configure_linesearch(DerivativeProbe(),nothing,x,Fx;value_derivative! = cb.joint!)
    @test isempty(cb.events)
    for mode in (:separate,:joint,:both)
        cb=observed_callbacks(T)
        search=configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;
            derivative! = mode===:joint ? nothing : cb.derivative!,
            value_derivative! = mode===:separate ? nothing : cb.joint!)
        r=perform_linesearch!(search,x,d,Fx)
        expected=mode===:separate ? (1,2) : mode===:joint ? (2,2) : (1,2)
        @test r.success && counts(r)==expected
        @test r.nfev==cb.calls[:values]+cb.calls[:joint]
        @test r.ndev==cb.calls[:derivative]+cb.calls[:joint]
        expected_events=mode===:separate ? [:derivative,:values,:derivative] :
            mode===:joint ? [:joint,:joint] : [:derivative,:joint]
        @test cb.events==expected_events
    end
    # Actual derivative-based rejections/reductions, with the same numerical result.
    for mode in (:separate,:joint,:both)
        cb=observed_callbacks(T;model=:quadratic)
        search=configure_linesearch(DerivativeProbe(),cb.value!,x,ones(T,2);
            derivative! = mode===:joint ? nothing : cb.derivative!,
            value_derivative! = mode===:separate ? nothing : cb.joint!)
        r=perform_linesearch!(search,x,T[3],ones(T,2))
        expected=mode===:separate ? (3,4) : mode===:joint ? (4,4) : (3,4)
        @test r.success && r.alpha==T(0.25) && r.iterations==2 && counts(r)==expected
        @test r.x_trial==T[0.75] && r.F_trial==fill(T(0.0625),2)
        @test counts(r)==observed_counts(cb)
    end
    # Errors and nonfinite outputs at base, value-only, derivative-only and joint boundaries.
    cases=((:separate,:derivative,1,:error,(0,1),:evaluation_failure),
           (:separate,:values,1,:error,(1,1),:evaluation_failure),
           (:separate,:derivative,2,:error,(1,2),:evaluation_failure),
           (:joint,:joint,1,:error,(1,1),:evaluation_failure),
           (:both,:joint,1,:error,(1,2),:evaluation_failure),
           (:separate,:values,1,:nonfinite_values,(1,1),:nonfinite_trial),
           (:separate,:derivative,1,:nonfinite_derivatives,(0,1),:nonfinite_trial),
           (:separate,:derivative,2,:nonfinite_derivatives,(1,2),:nonfinite_trial),
           (:joint,:joint,1,:nonfinite_values,(1,1),:nonfinite_trial),
           (:both,:joint,1,:nonfinite_values,(1,2),:nonfinite_trial),
           (:both,:joint,1,:nonfinite_derivatives,(1,2),:nonfinite_trial))
    for (mode,operation,call,action,expected,reason) in cases
        cb=observed_callbacks(T); cb.fault[]=(;operation,call,action)
        search=configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;
            derivative! = mode===:joint ? nothing : cb.derivative!,
            value_derivative! = mode===:separate ? nothing : cb.joint!)
        r=perform_linesearch!(search,x,d,Fx)
        @test !r.success && r.reason===reason && counts(r)==expected
        @test counts(r)==observed_counts(cb)
        reason===:evaluation_failure && @test r.message=="reported failure"
        @test x==zeros(T,1) && Fx==zeros(T,2)
    end
    for (operation,call,exception) in ((:values,1,ArgumentError),(:derivative,1,ErrorException),(:joint,1,ErrorException))
        cb=observed_callbacks(T)
        cb.fault[]=(;operation,call,action=exception===ArgumentError ? :invalid : :throw)
        search=configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;
            derivative! = cb.derivative!,value_derivative! = operation===:joint ? cb.joint! : nothing)
        @test_throws exception perform_linesearch!(search,x,d,Fx)
    end
    # Every component matters; declaring a capability is not a scalarization.
    positive_component!(D,y,d)=(fill!(D,-one(T));D[end]=one(T);nothing)
    search=configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;derivative! = positive_component!)
    r=perform_linesearch!(search,x,d,Fx)
    @test !r.success && r.reason===:not_descent && counts(r)==(0,1)
    # Container applicability is checked again when the caller changes direction type.
    vector_only!(D,y,d::Vector)=(fill!(D,-one(T));nothing)
    cb=observed_callbacks(T)
    search=configure_linesearch(DerivativeProbe(),cb.value!,x,Fx;derivative! = vector_only!)
    @test_throws ArgumentError perform_linesearch!(search,x,view(d,:),Fx)
    @test isempty(cb.events)
end

history_snapshot(s)=(copy(s.last_x),copy(s.last_F),copy(s.history[:,1:s.used]),s.head,s.used,s.accepted)

@testset "History commits and isolation ($T)" for T in (Float32,Float64)
    cb=observed_callbacks(T;model=:quadratic)
    x,Fx=zeros(T,1),ones(T,2)
    spec=HistoryProbe(capacity=3)
    stats=PreparationCounts()
    search=configure_linesearch(spec,cb.value!,x,Fx;counts=stats)
    storage=(search.x_trial,search.F_trial,search.last_x,search.last_F,search.history)
    @test search.accepted==0 && search.used==1 && search.history[:,1]==Fx
    @test search.last_x !== x && search.last_F !== Fx
    r=perform_linesearch!(search,x,T[3],Fx)
    @test r.success && r.nfev==2 && r.iterations==1
    @test search.accepted==1 && search.used==2 # rejected trial at 3 was not inserted
    copy_accept!(x,Fx,r)
    snapshot=history_snapshot(search)
    for action in (:error,:throw,:invalid,:nonfinite_values)
        cb.fault[]=(operation=:values,call=cb.calls[:values]+1,action=action)
        if action in (:throw,:invalid)
            @test_throws (action===:throw ? ErrorException : ArgumentError) perform_linesearch!(search,x,T[-0.25],Fx)
        else
            failed=perform_linesearch!(search,x,T[-0.25],Fx)
            @test !failed.success && failed.reason===(action===:error ? :evaluation_failure : :nonfinite_trial)
        end
        @test isequal(history_snapshot(search),snapshot)
        @test search.x_trial != search.last_x # scratch changed, committed base did not
    end
    cb.fault[]=(operation=:none,call=0,action=:ok)
    failed=perform_linesearch!(search,x,T[1],Fx;maxstep=0)
    @test !failed.success && failed.reason===:step_limit && failed.nfev==0
    @test isequal(history_snapshot(search),snapshot)
    # Finite exhaustion also leaves history uncommitted.
    failed=perform_linesearch!(search,x,T[1e6],Fx)
    @test !failed.success && failed.reason===:probe_iterations
    @test isequal(history_snapshot(search),snapshot)
    r2=perform_linesearch!(search,x,T[-0.25],Fx)
    @test r2.success && search.accepted==2 && search.used==3
    @test search.history[1,:]==T[1,0.25,0.0625] # no duplicate base insertion
    copy_accept!(x,Fx,r2)
    @test all(a===b for (a,b) in zip(storage,
        (search.x_trial,search.F_trial,search.last_x,search.last_F,search.history)))
    @test (stats.preparations,stats.conversions,stats.validations)==(1,1,1)
    before=length(cb.events)
    @test_throws ArgumentError perform_linesearch!(search,x .+ one(T),T[-0.25],Fx)
    @test_throws ArgumentError perform_linesearch!(search,x,T[-0.25],Fx .+ one(T))
    @test length(cb.events)==before
    # The stored window affects acceptance, including a value increase.
    cb2=observed_callbacks(T;model=:quadratic)
    a,b=zeros(T,1),ones(T,2)
    window=configure_linesearch(HistoryProbe(capacity=2),cb2.value!,a,b)
    rr=perform_linesearch!(window,a,T[1],b);copy_accept!(a,b,rr)
    rr=perform_linesearch!(window,a,T[0.75],b)
    @test rr.success && all(rr.F_trial .> b) && window.accepted==2
    @test window.used==2 && window.head==1 && window.history[1,:]==T[0.5625,0]
    # Interleaved preparations from the very same specification have disjoint state.
    va=observed_callbacks(T;model=:quadratic);vb=observed_callbacks(T;model=:quadratic)
    xa,fa=zeros(T,1),ones(T,2);xb,fb=T[4],fill(T(9),2)
    sa=configure_linesearch(spec,va.value!,xa,fa);sb=configure_linesearch(spec,vb.value!,xb,fb)
    @test all(!Base.mightalias(a,b) for a in ProtocolFixtures.outputs(sa), b in ProtocolFixtures.outputs(sb))
    ra=only(consume!(sa,xa,fa,[(T[3],(;))])); saved=history_snapshot(sa)
    point=copy(ra.x_trial);values=copy(ra.F_trial)
    rb=only(consume!(sb,xb,fb,[(T[-1],(;))]))
    @test isequal(history_snapshot(sa),saved) && ra.x_trial==point && ra.F_trial==values
    @test sa.accepted==sb.accepted==1 && xa!=xb && fa!=fb
    ra2=only(consume!(sa,xa,fa,[(T[-0.25],(;))]))
    @test ra2.success && sa.accepted==2 && sb.accepted==1
    @test rb.x_trial==xb && rb.F_trial==fb
end

# Source option conversion is observable without changing production code.
@testset "Preparation once and reused buffers ($T)" for T in (Float32,Float64)
    converted=Ref(0)
    wrap(v)=CountedReal(v,converted)
    spec=BackTracking(initial_step=InitialStatic(alpha=wrap(1.0)),c_1=wrap(1e-4),
        tau_1=wrap(0.05),tau_2=wrap(0.95),minstep=wrap(1e-15),maxstep=wrap(1.0))
    cb=observed_callbacks(T); x,Fx=zeros(T,1),zeros(T,2)
    search=configure_linesearch(spec,cb.value!,x,Fx)
    @test converted[]==6 && isempty(cb.events)
    buffers=(search.x_trial,search.F_trial)
    steps=[(T[0.25],(slopes=fill(T(-0.25),2),factor=T(-0.25))) for _ in 1:4]
    results=consume!(search,x,Fx,steps)
    @test all(r -> r.success,results) && converted[]==6
    @test all(r -> r.x_trial===buffers[1] && r.F_trial===buffers[2],results)
    for spec in (DerivativeProbe(),HistoryProbe())
        stats=PreparationCounts(); cb=observed_callbacks(T)
        x,Fx=zeros(T,1),zeros(T,2)
        search=configure_linesearch(spec,cb.value!,x,Fx;derivative! = cb.derivative!,counts=stats)
        buffers=ProtocolFixtures.outputs(search)
        results=consume!(search,x,Fx,[(T[0.25],(;)) for _ in 1:4])
        @test all(r -> r.success,results)
        @test (stats.preparations,stats.conversions,stats.validations)==(1,1,1)
        @test all(a===b for (a,b) in zip(buffers,ProtocolFixtures.outputs(search)))
        before=length(cb.events)
        @test_throws ArgumentError perform_linesearch!(search,x,T[NaN],Fx)
        @test length(cb.events)==before
    end
end

@testset "Interleaved executions of all methods ($T)" for T in (Float32,Float64)
    for spec in (BackTracking(),DerivativeProbe(),HistoryProbe())
        a=observed_callbacks(T);b=observed_callbacks(T)
        xa,fa=zeros(T,2),zeros(T,2);xb,fb=fill(T(3),2),fill(T(-6),2)
        sa=configure_linesearch(spec,a.value!,xa,fa;derivative! = a.derivative!)
        sb=configure_linesearch(spec,b.value!,xb,fb;derivative! = b.derivative!)
        owned_a=spec isa BackTracking ? (sa.x_trial,sa.F_trial) : ProtocolFixtures.outputs(sa)
        owned_b=spec isa BackTracking ? (sb.x_trial,sb.F_trial) : ProtocolFixtures.outputs(sb)
        @test all(!Base.mightalias(u,v) for u in owned_a,v in owned_b)
        data=spec isa BackTracking ? (slopes=-ones(T,2),factor=-one(T)) : (;)
        ra=only(consume!(sa,xa,fa,[(ones(T,2),data)]))
        saved=(copy(ra.x_trial),copy(ra.F_trial),scalars(ra))
        derivative_snapshot=spec isa DerivativeProbe ? (copy(sa.D_base),copy(sa.D_trial)) : nothing
        rb=only(consume!(sb,xb,fb,[(fill(T(0.5),2),data)]))
        @test ra.success && rb.success
        @test (ra.x_trial,ra.F_trial,scalars(ra))==saved
        if spec isa DerivativeProbe
            @test (sa.D_base,sa.D_trial)==derivative_snapshot
        end
        ra2=only(consume!(sa,xa,fa,[(fill(T(0.25),2),data)]))
        @test ra2.success && rb.x_trial==xb && rb.F_trial==fb
        if spec isa HistoryProbe
            @test sa.accepted==2 && sb.accepted==1
        end
    end
end

struct ShiftedVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(v::ShiftedVector)=size(v.data)
Base.axes(v::ShiftedVector)=(0:length(v.data)-1,)
Base.getindex(v::ShiftedVector,i::Int)=v.data[i+1]
Base.similar(v::ShiftedVector)=similar(v.data)

struct AliasedSimilar{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(v::AliasedSimilar)=size(v.data)
Base.getindex(v::AliasedSimilar,i::Int)=v.data[i]
Base.similar(v::AliasedSimilar)=v

@testset "Runtime structure, ownership and domain ($T)" for T in (Float32,Float64)
    for spec in (BackTracking(),DerivativeProbe(),HistoryProbe())
        cb=observed_callbacks(T);x,Fx=zeros(T,2),zeros(T,2);d=ones(T,2)
        data=spec isa BackTracking ? (slopes=-ones(T,2),factor=-one(T)) : (;)
        search=configure_linesearch(spec,cb.value!,x,Fx;derivative! = cb.derivative!)
        seed=perform_linesearch!(search,x,d,Fx;maxstep=0.25,data...)
        copy_accept!(x,Fx,seed)
        empty!(cb.events);empty!(cb.points)
        for key in keys(cb.calls)
            cb.calls[key]=0
        end
        for (a,b,c,error) in ((T[],T[],Fx,ArgumentError),(x,ones(T,3),Fx,DimensionMismatch),
                (x,d,zeros(T,3),DimensionMismatch),(zeros(Int,2),d,Fx,ArgumentError),
                (x,d,zeros(Int,2),ArgumentError),(fill(T(NaN),2),d,Fx,ArgumentError),
                (x,fill(T(Inf),2),Fx,ArgumentError),(x,d,fill(T(NaN),2),ArgumentError),
                (ShiftedVector(copy(x)),d,Fx,DimensionMismatch))
            @test_throws error perform_linesearch!(search,a,b,c;data...)
            @test isempty(cb.events)
        end
        for input in (search.x_trial,view(search.x_trial,:),search.F_trial,view(search.F_trial,:))
            @test_throws ArgumentError perform_linesearch!(search,input,d,Fx;data...)
            @test_throws ArgumentError perform_linesearch!(search,x,input,Fx;data...)
            @test_throws ArgumentError perform_linesearch!(search,x,d,input;data...)
        end
        @test isempty(cb.events)
        saved=copy(x);resize!(x,3)
        @test_throws DimensionMismatch perform_linesearch!(search,x,ones(T,3),Fx;data...)
        resize!(x,2);copyto!(x,saved)
        @test_throws DimensionMismatch configure_linesearch(spec,cb.value!,ShiftedVector(copy(x)),Fx;
            derivative! = cb.derivative!)
        @test_throws ArgumentError configure_linesearch(spec,cb.value!,AliasedSimilar(copy(x)),Fx;
            derivative! = cb.derivative!)
        for bound in (-1,NaN,"invalid")
            @test_throws ArgumentError perform_linesearch!(search,x,d,Fx;maxstep=bound,data...)
        end
        r=perform_linesearch!(search,x,d,Fx;maxstep=0,data...)
        @test !r.success && r.reason===:step_limit && counts(r)==(0,0) && isempty(cb.events)
        bound=BigFloat(T(0.5))-BigFloat(eps(T(0.5)))/4
        r=perform_linesearch!(search,view(x,:),view(d,:),view(Fx,:);maxstep=bound,data...)
        @test r.success && r.alpha===prevfloat(T(0.5)) && r.alpha<=bound
        @test all(y -> y==x || y==r.x_trial,cb.points)
    end
    for spec in (DerivativeProbe(alpha=1e-100),HistoryProbe(alpha=1e-100))
        if T===Float32
            @test_throws ArgumentError configure_linesearch(spec,(F,y)->nothing,T[0],T[0];
                derivative! = (D,y,d)->nothing)
        end
    end
    @test_throws ArgumentError configure_linesearch(HistoryProbe(capacity=0),(F,y)->nothing,T[0],T[0])
    for spec in (BackTracking(minstep=0.25),DerivativeProbe(minstep=0.25),HistoryProbe(minstep=0.25))
        cb=observed_callbacks(T);x,Fx=zeros(T,1),zeros(T,2)
        data=spec isa BackTracking ? (slopes=-ones(T,2),factor=-one(T)) : (;)
        search=configure_linesearch(spec,cb.value!,x,Fx;derivative! = cb.derivative!)
        r=perform_linesearch!(search,x,ones(T,1),Fx;maxstep=0.125,data...)
        @test !r.success && r.reason===:step_limit && r.alpha==T(0.125)
        @test counts(r)==(0,0) && isempty(cb.events)
    end
end
