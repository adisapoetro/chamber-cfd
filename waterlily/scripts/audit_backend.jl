# Independent engine-identity and analytical 3D checks; no chamber validation.
using ChamberWaterLily, WaterLily, Pkg, TOML, SHA, JSON3, LinearAlgebra
const CW = ChamberWaterLily

function abc_flow(n)
    # Periodic ABC flow is Beltrami: curl(u)=u and Laplacian(u)=-u on [0,2pi]^3.
    # Its exact viscous solution is u(t)=exp(-nu*t)u(0), p=-|u|^2/2.
    dx = 2pi/n
    viscosity = 0.05
    initial(i, x) = i == 1 ? sin(dx*x[3]) + cos(dx*x[2]) :
                    i == 2 ? sin(dx*x[1]) + cos(dx*x[3]) :
                             sin(dx*x[2]) + cos(dx*x[1])
    sim = WaterLily.Simulation((n,n,n), (0.,0.,0.), n; U=1.,
        ν=viscosity/dx, Δt=0.005/dx, u0=initial, perdir=(1,2,3), T=Float64,
        pois_ctor=f -> CW.CheckedPoisson(
            WaterLily.MultiLevelPoisson(f.p,f.μ₀,f.σ;perdir=(1,2,3)),1e-14,1000))
    for _ in 1:20
        sim.flow.Δt[end] = 0.005/dx
        WaterLily.sim_step!(sim;remeasure=false)
    end
    physical_time = WaterLily.time(sim)*dx
    error2 = 0.; exact2 = 0.; divergence = 0.
    for I in WaterLily.inside(sim.flow.p)
        divergence = max(divergence,abs(WaterLily.div(I,sim.flow.u))/dx)
        for k in 1:3
            exact = initial(k,WaterLily.loc(k,I,Float64))*exp(-viscosity*physical_time)
            error2 += (sim.flow.u[I,k]-exact)^2
            exact2 += exact^2
        end
    end
    return Dict("cells_per_axis"=>n,"physical_time_s"=>physical_time,
        "relative_velocity_l2_error"=>sqrt(error2/exact2),
        "max_divergence_s_inv"=>divergence,
        "all_three_velocity_components_nonzero"=>all(
            maximum(abs,@view(sim.flow.u[:,:,:,k]))>0 for k in 1:3))
end

function periodic_scalar_mesh(n)
    dims=(n,n,n); cells=vec(collect(CartesianIndices(dims)))
    ids=reshape(collect(1:length(cells)),dims)
    left=Int[];right=Int[];faces=CartesianIndex{3}[];direction=Int[]
    for (a,I) in enumerate(cells),k in 1:3
        J=CartesianIndex(ntuple(d -> d==k ? mod1(I[d]+1,n) : I[d],3))
        push!(left,a);push!(right,ids[J]);push!(faces,J);push!(direction,k)
    end
    dx=1/n;nf=length(left)
    return CW.TransportMesh(dims,cells,ids,left,right,faces,direction,
        ones(nf),fill(dx^2,nf),fill(dx^3,length(cells)),fill(dx^3,length(cells)),dx)
end

function scalar_wave(n)
    # Exact cell means of a Fourier mode transported in all three directions.
    # c = 400 + sinc(1/n)^3 exp(-3*(2pi)^2*D*t)
    #                 * sin(2pi*(x+y+z-(ux+uy+uz)*t)).
    mesh=periodic_scalar_mesh(n);velocity=(0.2,-0.1,0.15);D=0.02
    phase=[sum((Tuple(I).-0.5)./n) for I in mesh.cells]
    c=400 .+ sinc(1/n)^3 .* sin.(2pi.*phase)
    q=mesh.area.*[velocity[k] for k in mesh.direction]
    steps=2n;dt=0.1/steps;worst_budget=0.
    for _ in 1:steps
        c,a=CW.transport_step(mesh,c,q,fill(D,length(c)),zeros(length(c)),dt,400.)
        worst_budget=max(worst_budget,abs(a.budget_error_ppm_m3))
    end
    exact=400 .+ sinc(1/n)^3*exp(-3*(2pi)^2*D*0.1) .*
        sin.(2pi.*(phase.-sum(velocity)*0.1))
    return Dict("cells_per_axis"=>n,"physical_time_s"=>0.1,
        "rms_concentration_error_ppm"=>sqrt(sum(abs2,c.-exact)/length(c)),
        "mass_budget_max_ppm_m3"=>worst_budget,
        "minimum_ppm"=>minimum(c),"maximum_ppm"=>maximum(c))
end

function audit_backend(output)
    ispath(output) && error("Choose a fresh audit output")
    root=normpath(joinpath(@__DIR__,".."))
    manifest=TOML.parsefile(joinpath(root,"Manifest.toml"))
    declared=only(manifest["deps"]["WaterLily"])
    native_root=pkgdir(WaterLily)
    native_hash=bytes2hex(Pkg.GitTools.tree_hash(native_root))
    c=load_config(joinpath(root,"configs","current.toml"));state=build_flow(c)
    identity=Dict(
        "active_project_matches"=>normpath(Base.active_project())==joinpath(root,"Project.toml"),
        "julia_version"=>string(VERSION),"waterlily_version"=>string(pkgversion(WaterLily)),
        "declared_waterlily_tree_sha1"=>declared["git-tree-sha1"],
        "loaded_waterlily_tree_sha1"=>native_hash,
        "loaded_tree_matches_manifest"=>native_hash==declared["git-tree-sha1"],
        "backend"=>string(WaterLily.backend),
        "flow_type_is_native"=>state.sim.flow isa WaterLily.Flow,
        "pressure_type_is_native"=>state.sim.pois.inner isa WaterLily.MultiLevelPoisson,
        "scalar_transport"=>"Project-owned Julia finite-volume extension",
        "numerical_wall_m"=>CW.wall_thickness(c),
        "physical_wall_input_m"=>c["chamber"]["physical_wall_thickness_m"],
        "actual_domain_m"=>collect(state.dims).*state.sim.body.dx)
    airflow=[abc_flow(n) for n in (8,16,32)]
    scalar=[scalar_wave(n) for n in (8,12,16)]
    ue=[x["relative_velocity_l2_error"] for x in airflow]
    ce=[x["rms_concentration_error_ppm"] for x in scalar]
    checks=Dict(
        "loaded_registered_package"=>identity["loaded_tree_matches_manifest"] &&
            identity["waterlily_version"]=="1.8.0" && identity["active_project_matches"],
        "native_flow_and_pressure"=>identity["flow_type_is_native"] && identity["pressure_type_is_native"],
        "three_dimensional_airflow_refines"=>all(diff(ue).<0) && ue[end]<0.005,
        "three_dimensional_airflow_divergence"=>maximum(x["max_divergence_s_inv"] for x in airflow)<1e-5,
        "three_dimensional_transport_refines"=>all(diff(ce).<0) && ce[end]<0.1,
        "three_dimensional_transport_conserves"=>maximum(x["mass_budget_max_ppm_m3"] for x in scalar)<1e-9,
        "three_dimensional_transport_bounded"=>all(399<=x["minimum_ppm"]<=x["maximum_ppm"]<=401 for x in scalar))
    sources=Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(root,"src")),name in files
        path=joinpath(dir,name);sources[relpath(path,root)]=bytes2hex(sha256(read(path)))
    end
    for name in ["scripts/audit_backend.jl","Project.toml","Manifest.toml","LocalPreferences.toml","configs/current.toml"]
        sources[name]=bytes2hex(sha256(read(joinpath(root,name))))
    end
    report=Dict("status"=>all(values(checks)) ? "PASS_IMPLEMENTATION_BENCHMARKS" : "FAIL_IMPLEMENTATION_BENCHMARKS",
        "scientific_scope"=>"Canonical periodic flows and scalar transport only; no chamber, moving-wall, fan or field validation",
        "identity"=>identity,"native_3d_abc_flow"=>airflow,"extension_3d_advection_diffusion"=>scalar,
        "checks"=>checks,"source_sha256"=>sources)
    mkpath(dirname(output));CW.write_json(output,report)
    println(JSON3.write(report))
    return all(values(checks))
end

length(ARGS)==1 || error("Usage: julia --project=. scripts/audit_backend.jl results/engine_audit.json")
audit_backend(ARGS[1]) || exit(1)
