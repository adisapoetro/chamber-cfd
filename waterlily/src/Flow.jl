"""Delegate pressure solves to WaterLily with explicit convergence controls."""
struct CheckedPoisson{T,A,V,P<:WaterLily.AbstractPoisson{T,A,V}} <: WaterLily.AbstractPoisson{T,A,V}
    inner::P
    tolerance::Float64
    iterations::Int
end
CheckedPoisson(p::WaterLily.AbstractPoisson{T,A,V},tol,it) where {T,A,V}=CheckedPoisson{T,A,V,typeof(p)}(p,tol,it)
Base.getproperty(p::CheckedPoisson,k::Symbol)=k in (:inner,:tolerance,:iterations) ? getfield(p,k) : getproperty(getfield(p,:inner),k)
WaterLily.update!(p::CheckedPoisson)=WaterLily.update!(p.inner)
function WaterLily.solver!(p::CheckedPoisson;kwargs...)
    WaterLily.solver!(p.inner;tol=p.tolerance,itmx=p.iterations,kwargs...)
    r=WaterLily.L₂(p.inner.levels[1])
    isfinite(r) && r<=p.tolerance || error("Native WaterLily pressure solve failed: residual²=$r")
end

function build_flow(c)
    validate_config(c)
    body=chamber_body(c);dx=body.dx;nu=c["numerics"]["viscosity_m2_s"]
    dims=Tuple(2ceil(Int,c["domain"][k]/(2dx)) for k in ["x_m","y_m","z_m"])
    speed=c["air"]["wind_speed_m_s"];angle=c["air"]["wind_to_deg"]
    wind=SVector(speed*cosd(angle),speed*sind(angle),0.)
    bc(i,x,t)=body.lower[3]+x[3]*dx<0 ? 0. : wind[i]
    function initial(i,x)
        p=body.lower+x*dx;ch=c["chamber"]
        inside=abs(p[1])<=ch["width_m"]/2 && abs(p[2])<=ch["depth_m"]/2 && ch["floor_m"]<=p[3]<=ch["floor_m"]+ch["height_m"]
        return inside || WaterLily.sdf(body,x,0.)<=0 ? 0. : wind[i]
    end
    w=c["waterlily"]
    sim=WaterLily.Simulation(dims,bc,c["chamber"]["width_m"]/dx;U=1.,ν=nu/dx,
        Δt=min(c["numerics"]["dt_s"]/dx,0.1),ϵ=w["kernel_width_cells"],body,u0=initial,T=Float64,exitBC=false,
        pois_ctor=flow->CheckedPoisson(WaterLily.MultiLevelPoisson(flow.p,flow.μ₀,flow.σ),w["pressure_tolerance"],w["pressure_max_iterations"]))
    forces=forcing_fields(c,sim)
    volumes=reference_volumes(c,dims,body.lower)
    isapprox(sum(volumes),nominal_volume(c);rtol=1e-12) || error("Reference-volume mismatch")
    measured_motion=Ref((0.,0.))
    return (;sim,c,dims,forces,volumes,measured_motion)
end

function advance_flow!(state,target_s;after_step=nothing)
    sim=state.sim;c=state.c;dx=sim.body.dx;f=state.forces
    while WaterLily.time(sim)*dx<target_s-1e-9
        t=WaterLily.time(sim)*dx
        rate_dt=f.max_rate[]>0 ? 0.4/f.max_rate[] : Inf
        visc_dt=0.4/(6(sim.flow.ν+f.max_eddy[])+eps())
        dt=min(sim.flow.Δt[end],c["numerics"]["dt_s"]/dx,rate_dt,visc_dt,(target_s-t)/dx,(next_event(c,t)-t)/dx)
        dt>1e-12 && isfinite(dt) || error("Invalid time step")
        sim.flow.Δt[end]=dt
        f.active[]=cycle_state(c,t+dt*dx/2).fans_on
        # The native predictor/corrector advances momentum and pressure.
        next_motion=motion(sim.body,t+dt*dx)
        remeasure=next_motion!=state.measured_motion[]
        WaterLily.sim_step!(sim;remeasure,udf=f.force!)
        state.measured_motion[]=next_motion
        all(isfinite,sim.flow.u) || error("Non-finite WaterLily velocity")
        isnothing(after_step) || after_step(state,t,t+dt*dx)
    end
    return state
end

function flow_observables(state)
    sim=state.sim;dx=sim.body.dx;t=WaterLily.time(sim)*dx
    maxdiv=0.;wallerror=0.;corediv=0.;speed=zeros(state.dims)
    for I in WaterLily.inside(sim.flow.p)
        d,n,v=WaterLily.measure(sim.body,WaterLily.loc(0,I,Float64),t/dx)
        velocity=SVector(ntuple(k->(sim.flow.u[I,k]+sim.flow.u[I+WaterLily.δ(k,I),k])/2,3))
        div=abs(WaterLily.div(I,sim.flow.u))/dx
        d>sim.ϵ && (maxdiv=max(maxdiv,div))
        d < -sim.ϵ && (wallerror=max(wallerror,norm(velocity-v));corediv=max(corediv,div))
        speed[I-CartesianIndex(1,1,1)]=norm(velocity)
    end
    return Dict{String,Any}("time_s"=>t,"gap_m"=>cycle_state(state.c,t).gap,"fans_on"=>Int(cycle_state(state.c,t).fans_on),
        "roi_volume_m3"=>sum(state.volumes),"roi_mean_speed_m_s"=>sum(speed.*state.volumes)/sum(state.volumes),
        "max_speed_m_s"=>maximum(speed),"fluid_divergence_max_s_inv"=>maxdiv,"solid_core_velocity_error_m_s"=>wallerror,
        "solid_core_divergence_max_s_inv"=>corediv,"pressure_residual_squared"=>WaterLily.L₂(sim.pois.inner.levels[1]))
end

function save_fields(state,out,index)
    stem=@sprintf("%04d",index);sim=state.sim
    record=flow_observables(state)
    record["velocity_file"]="fields/$(stem)_velocity.f64"
    record["distance_file"]="fields/$(stem)_distance.f64"
    vel=zeros(Float64,state.dims...,3);dist=zeros(Float64,state.dims)
    for I in CartesianIndices(dist)
        J=I+CartesianIndex(1,1,1)
        for k in 1:3
            vel[I,k]=(sim.flow.u[J,k]+sim.flow.u[J+WaterLily.δ(k,J),k])/2
        end
        dist[I]=sim.body.dx*WaterLily.sdf(sim.body,WaterLily.loc(0,J,Float64),record["time_s"]/sim.body.dx)
    end
    open(io->write(io,vel),joinpath(out,record["velocity_file"]),"w")
    open(io->write(io,dist),joinpath(out,record["distance_file"]),"w")
    return record
end

function write_json(path,value)
    open(path,"w") do io
        JSON3.pretty(io,value);println(io)
    end
end

function run_airflow(c,out;duration_s=c["cycle"]["cycles"]*period(c))
    0<duration_s<=c["cycle"]["cycles"]*period(c) || error("Duration is outside the configured cycle record")
    ispath(out) && error("Choose a fresh output folder: $out")
    mkpath(joinpath(out,"fields"));state=build_flow(c)
    write_json(joinpath(out,"inputs.json"),c)
    root=normpath(joinpath(@__DIR__,".."))
    source=Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(root,"src")),name in files
        path=joinpath(dir,name);source[relpath(path,root)]=bytes2hex(sha256(read(path)))
    end
    for name in ["Project.toml","Manifest.toml","LocalPreferences.toml","scripts/run.jl","scripts/render_airflow.py",
                 "configs/current.toml","configs/reference_python_input.yaml","assets/reference_palm.json","requirements-render.lock"]
        source[name]=bytes2hex(sha256(read(joinpath(root,name))))
    end
    for path in keys(source)
        destination=joinpath(out,"configs","source",path)
        mkpath(dirname(destination));cp(joinpath(root,path),destination)
    end
    write_json(joinpath(out,"configs","source_hashes.json"),source)
    write_json(joinpath(out,"progress.json"),Dict("status"=>"RUNNING","time_s"=>0.,"target_s"=>duration_s,"co2_transport_enabled"=>false))
    records=Any[];push!(records,save_fields(state,out,0));start=time()
    saves=collect(c["numerics"]["save_every_s"]:c["numerics"]["save_every_s"]:duration_s)
    (isempty(saves)||saves[end]<duration_s) && push!(saves,duration_s)
    for (i,t) in enumerate(saves)
        advance_flow!(state,t);record=save_fields(state,out,i);push!(records,record)
        println(JSON3.write(record));flush(stdout)
        write_json(joinpath(out,"progress.json"),Dict("status"=>"RUNNING","time_s"=>t,"target_s"=>duration_s,"co2_transport_enabled"=>false))
    end
    result=Dict("engine"=>"WaterLily","version"=>string(pkgversion(WaterLily)),"julia"=>string(VERSION),
        "scientific_status"=>"DIAGNOSTIC_NOT_VALIDATED","co2_transport_enabled"=>false,
        "waterlily_backend"=>WaterLily.backend,"julia_threads"=>Threads.nthreads(),"blas_threads"=>LinearAlgebra.BLAS.get_num_threads(),
        "exchange_input_applied"=>false,"elapsed_s"=>time()-start,"duration_s"=>duration_s,
        "source_sha256"=>source,"shape"=>state.dims,"array_order"=>"Fortran/Julia column-major",
        "grid_cell_m"=>state.sim.body.dx,"grid_lower_m"=>state.sim.body.lower,
        "actual_domain_m"=>collect(state.dims).*state.sim.body.dx,
        "nominal_reference_volume_m3"=>nominal_volume(c),"numerical_wall_thickness_m"=>wall_thickness(c),
        "reference_input_exchange_umol_s"=>c["exchange"]["net_co2_umol_s"],
        "fan_configuration"=>fan_configuration(c),"snapshots"=>records,
        "limitations"=>["Airflow only; no CO2 field or mixing-time inference.","Outward wall thickening changes the exterior geometry.",
                        "Native BDIM wall behaviour needs refinement checks.","Native external velocity boundary conditions differ from the Python reference.",
                        "Plant drag is schematic; fan performance and support position are assumed."])
    write_json(joinpath(out,"summary.json"),result)
    write_json(joinpath(out,"progress.json"),Dict("status"=>"COMPLETE_AIRFLOW_ONLY","time_s"=>duration_s,"co2_transport_enabled"=>false))
    names=["time_s","gap_m","fans_on","roi_volume_m3","roi_mean_speed_m_s","max_speed_m_s","fluid_divergence_max_s_inv","solid_core_velocity_error_m_s"]
    open(joinpath(out,"timeseries.csv"),"w") do io
        println(io,join(names,","))
        for r in records;println(io,join([r[k] for k in names],","));end
    end
    hashes=Dict{String,String}()
    for (dir,_,files) in walkdir(out),name in files
        p=joinpath(dir,name);hashes[relpath(p,out)]=bytes2hex(sha256(read(p)))
    end
    write_json(joinpath(out,"manifest.json"),Dict("sha256"=>hashes))
    return result
end
