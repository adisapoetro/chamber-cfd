mutable struct CoupledScalar
    concentration::Array{Float64,3}
    leaf::Any
    mesh::TransportMesh
    projection::Any
    gap::Float64
    audits::Vector{Any}
    boundary_inventory::Float64
    source_inventory::Float64
    initial_inventory::Float64
end

function build_scalar(state)
    b=state.sim.body;c=state.c
    mesh=transport_mesh(b,state.dims,0.,.01);leaf=leaflet_source(state)
    C=fill(c["air"]["ambient_co2_ppm"],state.dims)
    labels,exposed=mesh_components(mesh)
    center=SVector(0.,0.,b.floor+b.size[3]/2)
    seed=CartesianIndex(Tuple(floor.(Int,(center-b.lower)/b.dx).+1))
    chamber=labels[mesh.ids[seed]]
    chamber in exposed && error("Closed chamber is connected to ambient")
    interior=findall(==(chamber),labels)
    isapprox(sum(mesh.old_volume[interior]),nominal_volume(c);atol=1e-9) || error("Closed capacity differs from nominal chamber")
    C[mesh.cells[interior]].=c["air"]["initial_co2_ppm"]
    all(I in mesh.cells[interior] for I in findall(!iszero,leaf.source)) || error("Leaf source is outside the closed airspace")
    initial=sum(mesh.old_volume.*C[mesh.cells])
    return CoupledScalar(C,leaf,mesh,prepare_projection(mesh),0.,Any[],0.,0.,initial)
end

function advance_scalar!(scalar,state,t0,t1)
    dt=t1-t0;b=state.sim.body;c=state.c
    gap0,_=motion(b,t0);gap1,_=motion(b,t1)
    if !(gap0==gap1==scalar.gap)
        scalar.mesh=transport_mesh(b,state.dims,t0,t1)
        scalar.projection=prepare_projection(scalar.mesh)
        scalar.gap=gap0==gap1 ? gap1 : NaN
    end
    mesh=scalar.mesh
    raw,mapping=waterlily_face_flux(state,mesh,t0,t1)
    q,a=compatible_flux(mesh,raw,dt;prepared=scalar.projection)
    num=c["numerics"];Cs=Float64(num["smagorinsky"]);S=state.forces.S
    D=[num["molecular_diffusivity_m2_s"]+b.dx*smagorinsky(I+CartesianIndex(1,1,1);S,Cs,Δ=1.)/num["turbulent_schmidt"] for I in mesh.cells]
    source=scalar.leaf.source[mesh.cells]
    isapprox(sum(source),sum(scalar.leaf.source);atol=1e-12) || error("Moving wall covers a leaflet source")
    cold=scalar.concentration[mesh.cells]
    C,z=transport_step(mesh,cold,q,D,source,dt,c["air"]["ambient_co2_ppm"])
    scalar.concentration[mesh.cells]=C
    scalar.boundary_inventory+=dt*z.boundary_outward_ppm_m3_s
    scalar.source_inventory+=dt*z.source_ppm_m3_s
    accumulated=sum(mesh.new_volume.*C)-scalar.initial_inventory-scalar.source_inventory+scalar.boundary_inventory
    abs(accumulated)<1e-5 || error("Cumulative CO₂ budget failed: $accumulated ppm m³")
    roi=[state.volumes[mesh.cells[mesh.left[e]]]>0 || (mesh.right[e]>0 && state.volumes[mesh.cells[mesh.right[e]]]>0) for e in eachindex(q)]
    interface=[mesh.area[e]<b.dx^2*(1-1e-10) for e in eachindex(q)]
    ratio(mask)=norm((q.-raw)[mask])/max(norm(raw[mask]),1e-30)
    audit=merge((time_s=t1,dt_s=dt),a,z,mapping,(cumulative_budget_error_ppm_m3=accumulated,
        roi_correction_relative_l2=ratio(roi),interface_correction_relative_l2=ratio(interface),
        roi_raw_flux_l2_m3_s=norm(raw[roi]),roi_delta_flux_l2_m3_s=norm((q.-raw)[roi]),
        co2_min_ppm=minimum(C),co2_max_ppm=maximum(C),max_diffusivity_m2_s=maximum(D)))
    push!(scalar.audits,audit)
    return scalar
end

function save_coupled(state,scalar,out,index)
    r=save_fields(state,out,index);stem=@sprintf("%04d",index)
    r["co2_file"]="fields/$(stem)_co2.f64";r["volume_file"]="fields/$(stem)_fluid_volume.f64"
    volume=zeros(state.dims);volume[scalar.mesh.cells]=scalar.mesh.new_volume
    open(io->write(io,scalar.concentration),joinpath(out,r["co2_file"]),"w")
    open(io->write(io,volume),joinpath(out,r["volume_file"]),"w")
    C=scalar.concentration;w=state.volumes
    mean=sum(C.*w)/sum(w);sd=sqrt(sum(w.*(C.-mean).^2)/sum(w))
    r["co2_mean_ppm"]=mean;r["co2_spatial_sd_ppm"]=sd
    r["co2_min_ppm"]=minimum(C[volume.>0]);r["co2_max_ppm"]=maximum(C[volume.>0])
    r["net_exchange_umol_s"]=state.c["exchange"]["net_co2_umol_s"]
    if !isempty(scalar.audits)
        r["cumulative_budget_error_ppm_m3"]=scalar.audits[end].cumulative_budget_error_ppm_m3
        r["roi_flux_correction_relative_l2"]=scalar.audits[end].roi_correction_relative_l2
    end
    return r
end

function freeze_source(out)
    root=normpath(joinpath(@__DIR__,".."));paths=String[]
    for sub in ["src","configs","assets","scripts","test","docs"]
        for (dir,_,files) in walkdir(joinpath(root,sub)),file in files
            endswith(file,".pyc") && continue
            push!(paths,relpath(joinpath(dir,file),root))
        end
    end
    append!(paths,["Project.toml","Manifest.toml","LocalPreferences.toml","requirements-render.lock","README.md"])
    hashes=Dict{String,String}()
    for p in paths
        destination=joinpath(out,"configs","source",p);mkpath(dirname(destination))
        cp(joinpath(root,p),destination);hashes[p]=bytes2hex(sha256(read(destination)))
    end
    write_json(joinpath(out,"configs","source_hashes.json"),hashes)
    return hashes
end

function output_manifest(out)
    hashes=Dict{String,String}()
    for (dir,_,files) in walkdir(out),name in files
        name=="manifest.json" && continue
        p=joinpath(dir,name);hashes[relpath(p,out)]=bytes2hex(sha256(read(p)))
    end
    write_json(joinpath(out,"manifest.json"),Dict("sha256"=>hashes))
end

function run_coupled(c,out;duration_s=c["cycle"]["cycles"]*period(c))
    0<duration_s<=c["cycle"]["cycles"]*period(c) || error("Invalid duration")
    ispath(out) && error("Choose a fresh output folder: $out")
    mkpath(joinpath(out,"fields"));state=build_flow(c);scalar=build_scalar(state)
    source=freeze_source(out);write_json(joinpath(out,"inputs.json"),c)
    open(io->write(io,scalar.leaf.area),joinpath(out,"fields","leaflet_area.f64"),"w")
    records=Any[save_coupled(state,scalar,out,0)];start=time()
    saves=collect(c["numerics"]["save_every_s"]:c["numerics"]["save_every_s"]:duration_s)
    (isempty(saves)||saves[end]<duration_s) && push!(saves,duration_s)
    for (i,t) in enumerate(saves)
        advance_flow!(state,t;after_step=(state,a,z)->advance_scalar!(scalar,state,a,z))
        record=save_coupled(state,scalar,out,i);push!(records,record)
        println(JSON3.write(record));flush(stdout)
        write_json(joinpath(out,"progress.json"),Dict("status"=>"RUNNING","time_s"=>t,"target_s"=>duration_s,"co2_transport_enabled"=>true))
    end
    audits=scalar.audits
    result=Dict{String,Any}("engine"=>"WaterLily","version"=>string(pkgversion(WaterLily)),"julia"=>string(VERSION),
        "scientific_status"=>"DIAGNOSTIC_NOT_VALIDATED","co2_transport_enabled"=>true,"exchange_input_applied"=>true,
        "transport"=>"Separate Julia conservative cut-cell finite-volume extension; implicit first-order upwind advection and diffusion",
        "waterlily_backend"=>WaterLily.backend,"julia_threads"=>Threads.nthreads(),"blas_threads"=>LinearAlgebra.BLAS.get_num_threads(),
        "elapsed_s"=>time()-start,"duration_s"=>duration_s,"source_sha256"=>source,"shape"=>state.dims,
        "array_order"=>"Fortran/Julia column-major","grid_cell_m"=>state.sim.body.dx,"grid_lower_m"=>state.sim.body.lower,
        "nominal_reference_volume_m3"=>nominal_volume(c),"numerical_wall_thickness_m"=>wall_thickness(c),
        "fan_configuration"=>fan_configuration(c),"snapshots"=>records,
        "leaf_source"=>Dict("triangles"=>scalar.leaf.triangles,"representative_area_m2"=>scalar.leaf.total_area_m2,
            "air_molar_density_mol_m3"=>scalar.leaf.molar_density,"net_exchange_umol_s"=>c["exchange"]["net_co2_umol_s"]),
        "checks"=>Dict("gcl_max_m3_s"=>maximum(a.gcl_max_m3_s for a in audits),
            "cumulative_budget_error_max_ppm_m3"=>maximum(abs(a.cumulative_budget_error_ppm_m3) for a in audits),
            "co2_min_ppm"=>minimum(a.co2_min_ppm for a in audits),"co2_max_ppm"=>maximum(a.co2_max_ppm for a in audits),
            "flux_correction_max_relative_l2"=>maximum(a.correction_relative_l2 for a in audits),
            "roi_flux_correction_max_relative_l2"=>maximum(a.roi_correction_relative_l2 for a in audits),
            "aperture_velocity_correction_max_m_s"=>maximum(a.correction_max_m_s for a in audits),
            "one_sided_extension_faces_max"=>maximum(a.one_sided_extension_faces for a in audits),
            "one_sided_extension_area_fraction_max"=>maximum(a.one_sided_extension_area_fraction for a in audits),
            "coupling_correction_10pct_gate"=>all(a.roi_correction_relative_l2<=.1 for a in audits) ? "PASS" : "FAIL"),
        "limitations"=>["Engineering diagnostic; no experimental validation or grid-convergence claim.",
            "Numerically thickened WaterLily walls alter external flow geometry.",
            "Sharp scalar apertures use endpoint-paired BDIM deblending and an audited auxiliary flux projection.",
            "Covered or newly exposed scalar faces use a counted one-sided velocity extension; unresolved closing gaps remain a limitation.",
            "First-order upwind transport adds numerical diffusion; refinement is required for mixing claims.",
            "Fan thrust, canopy drag and support placement are prescribed approximations.",
            "February 2025 net exchange is prescribed through both phases; no dynamic photosynthesis model."])
    write_json(joinpath(out,"summary.json"),result)
    write_json(joinpath(out,"progress.json"),Dict("status"=>"COMPLETE_COUPLED_DIAGNOSTIC","time_s"=>duration_s,"co2_transport_enabled"=>true))
    open(joinpath(out,"transport_audit.csv"),"w") do io
        println(io,join(string.(keys(first(audits))),","))
        for a in audits;println(io,join(values(a),","));end
    end
    names=["time_s","gap_m","fans_on","co2_mean_ppm","co2_spatial_sd_ppm","co2_min_ppm","co2_max_ppm","net_exchange_umol_s","roi_mean_speed_m_s"]
    open(joinpath(out,"timeseries.csv"),"w") do io
        println(io,join(names,","))
        for r in records;println(io,join([r[k] for k in names],","));end
    end
    output_manifest(out)
    return result
end
