using Test, ChamberWaterLily, WaterLily, StaticArrays, LinearAlgebra, TOML, JSON3, SHA
const CW=ChamberWaterLily
const INPUT=joinpath(@__DIR__,"..","configs","current.toml")
include("transport.jl")

@testset "Reference inputs and fan configuration" begin
    c=load_config(INPUT)
    @test nominal_volume(c)==96.
    @test c["exchange"]["net_co2_umol_s"]==-44.40047791098563
    @test c["palm"]["height_m"]==2.39
    fans=fan_configuration(c)
    @test [f.height for f in fans]≈[2/3,2,10/3]
    @test [f.azimuth for f in fans]==[0,120,240]
    @test norm(sum(f.direction*f.thrust for f in fans))<1e-14
    @test fans[1].thrust≈1.5383767444688012/3
    @test fans[1].position≈SVector(.3,0.,.6+2/3)
    for n in [0,1,4]
        q=deepcopy(c);q["fans"]["count"]=n
        @test length(fan_configuration(validate_config(q)))==n
    end
    q=deepcopy(c);q["fans"]["count"]=-1;@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["fans"]["heights_m"]=[1.,2.];@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["cycle"]["open_phase_s"]=10.;@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["chamber"]["height_m"]=2.;@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["waterlily"]["numerical_wall_cells"]=1.;@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["air"]["wind_speed_m_s"]=NaN;@test_throws ErrorException validate_config(q)
    q=deepcopy(c);q["air"]["windspeed"]=1.;@test_throws ErrorException validate_config(q)
end

@testset "Self-contained airflow bundle" begin
    c=load_config(INPUT)
    mktempdir() do tmp
        out=joinpath(tmp,"run")
        q=deepcopy(c)
        q["cycle"]["closed_s"]=.1;q["cycle"]["open_phase_s"]=.3
        q["cycle"]["opening_travel_s"]=.1;q["cycle"]["closing_travel_s"]=.1
        q["numerics"]["save_every_s"]=.1
        result=run_airflow(q,out;duration_s=.4)
        @test result["co2_transport_enabled"]===false
        @test result["exchange_input_applied"]===false
        @test result["duration_s"]==.4
        @test length(result["snapshots"])==5
        @test result["blas_threads"]==1
        parsed=JSON3.read(read(joinpath(out,"summary.json"),String))
        @test parsed.snapshots[end].time_s≈.4 atol=1e-12
        for (path,hash) in JSON3.read(read(joinpath(out,"manifest.json"),String)).sha256
            @test bytes2hex(sha256(read(joinpath(out,String(path)))))==hash
        end
        @test_throws ErrorException run_airflow(c,out;duration_s=.05)
    end
end

@testset "Five-minute closure and ten-minute open phase" begin
    c=load_config(INPUT)
    for (t,g,phase,fan) in [(0.,0.,"sealed",true),(299.,0.,"sealed",true),
            (300.,0.,"opening",false),(310.,1.,"opening",false),(320.,2.,"open",false),
            (879.,2.,"open",false),(880.,2.,"closing",false),(890.,1.,"closing",false),(900.,0.,"sealed",true)]
        s=cycle_state(c,t)
        @test s.gap≈g
        @test s.phase==phase
        @test s.fans_on==fan
        @test cycle_state(c,t+900).gap≈g
    end
    @test CW.next_event(c,299.)==300.
    @test CW.next_event(c,300.)==320.
    @test CW.next_event(c,899.5)==900.
end

@testset "Native body distance, wall velocity and reference volume" begin
    c=load_config(INPUT);b=chamber_body(c)
    grid(x)=(SVector{3,Float64}(x)-b.lower)/b.dx
    @test WaterLily.sdf(b,grid([0.,0.,2.6]),0.)>0
    @test WaterLily.sdf(b,grid([0.,0.,4.6]),0.)≈0 atol=1e-12
    @test WaterLily.sdf(b,grid([0.,0.,5.1]),0.)<0
    @test WaterLily.sdf(b,grid([0.,0.,5.1]),320/b.dx)>0
    @test WaterLily.sdf(b,grid([0.,0.,-.1]),0.)<0
    for (t,sign,v) in [(310.,-1.,-.05),(310.,1.,.05),(890.,-1.,.05),(890.,1.,-.05)]
        g=cycle_state(c,t).gap
        d,n,velocity=WaterLily.measure(b,grid([sign*(2+g/2),0.,2.6]),t/b.dx)
        @test d≈0 atol=1e-12
        @test velocity≈SVector(v,0.,0.)
        @test n≈SVector(-sign,0.,0.)
    end
    for dx in [.5,.4,.25]
        q=deepcopy(c);q["numerics"]["cell_m"]=dx;body=chamber_body(q)
        dims=Tuple(ceil(Int,q["domain"][k]/dx) for k in ["x_m","y_m","z_m"])
        @test sum(CW.reference_volumes(q,dims,body.lower))≈96. atol=1e-11
    end
end

@testset "Native WaterLily canonical flow" begin
    @test pkgversion(WaterLily)==v"1.8.0"
    @test !isdefined(WaterLily,:transport!)
    # A periodic Taylor-Green vortex has an independent exponential decay law.
    errors=Float64[]
    for n in [16,32]
        dx=2pi/n;k=dx;nu=.05/dx
        initial(i,x)=i==1 ? sin(k*x[1])*cos(k*x[2]) : -cos(k*x[1])*sin(k*x[2])
        sim=Simulation((n,n),(0.,0.),n;U=1.,ν=nu,Δt=.005/dx,u0=initial,perdir=(1,2),T=Float64,
            pois_ctor=flow->CW.CheckedPoisson(WaterLily.MultiLevelPoisson(flow.p,flow.μ₀,flow.σ;perdir=(1,2)),1e-14,1000))
        for _ in 1:40
            sim.flow.Δt[end]=.005/dx;sim_step!(sim;remeasure=false)
        end
        t=WaterLily.time(sim)
        err=0.;denom=0.
        for I in WaterLily.inside(sim.flow.p),d in 1:2
            expected=initial(d,WaterLily.loc(d,I,Float64))*exp(-2nu*k^2*t)
            err+=(sim.flow.u[I,d]-expected)^2;denom+=expected^2
        end
        push!(errors,sqrt(err/denom))
        @test all(isfinite,sim.flow.u)
        @test maximum(abs(WaterLily.div(I,sim.flow.u)) for I in WaterLily.inside(sim.flow.p))<1e-4
    end
    @test errors[2]<errors[1]
    @test errors[2]<.005
end

@testset "Native chamber flow and fan forcing" begin
    c=load_config(INPUT);state=build_flow(c)
    @test state.sim.flow isa WaterLily.Flow
    @test state.sim.pois.inner isa WaterLily.MultiLevelPoisson
    @test sum(state.volumes)≈96.
    field=state.forces.fan_field;dx=c["numerics"]["cell_m"]
    rho=c["air"]["pressure_pa"]/(287.05*c["air"]["temperature_k"])
    resultant=SVector(ntuple(i->sum(@view field[:,:,:,i])*dx^2*rho,3))
    @test norm(resultant)<1e-13
    @test maximum(abs,field)>0
    advance_flow!(state,.1)
    @test WaterLily.time(state.sim)*dx≈.1 atol=1e-12
    @test all(isfinite,state.sim.flow.u)
    stats=CW.flow_observables(state)
    @test stats["roi_mean_speed_m_s"]>0
    @test stats["pressure_residual_squared"]<=c["waterlily"]["pressure_tolerance"]
    # A fan-only control checks the source without prescribed external wind.
    q=deepcopy(c);q["air"]["wind_speed_m_s"]=0.
    still=build_flow(q);advance_flow!(still,.1)
    @test CW.flow_observables(still)["roi_mean_speed_m_s"]>0
end
