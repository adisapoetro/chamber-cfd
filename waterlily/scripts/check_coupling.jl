# Compare old/new mappings on identical native airflow; this does not solve CO₂.
using ChamberWaterLily, WaterLily, LinearAlgebra, JSON3, SHA
const CW=ChamberWaterLily

function legacy_midpoint_flux(state,mesh,t0,t1)
    sim=state.sim;b=sim.body;raw=zeros(length(mesh.left))
    for e in eachindex(raw)
        I=mesh.faces[e]+CartesianIndex(1,1,1);k=mesh.direction[e]
        d,_,v=WaterLily.measure(b,WaterLily.loc(k,I,Float64),(t0+t1)/(2b.dx))
        mu=WaterLily.μ₀(d,sim.ϵ)
        mu>1e-8 || error("Historical mapping lacks native support")
        raw[e]=mesh.orientation[e]*mesh.area[e]*(sim.flow.u[I,k]-(1-mu)*v[k])/mu
    end
    return raw
end

length(ARGS)==1 || error("Give a fresh JSON output path")
ispath(ARGS[1]) && error("Output exists")
c=load_config(joinpath(@__DIR__,"..","configs","current.toml"))
# Preserve real travel speed, but shorten the stationary periods for this probe.
c["cycle"]["closed_s"]=1.;c["cycle"]["open_phase_s"]=41.
c["cycle"]["opening_travel_s"]=20.;c["cycle"]["closing_travel_s"]=20.
state=build_flow(c);records=Any[]
elapsed=@elapsed advance_flow!(state,42.;after_step=function(s,t0,t1)
    CW.motion(s.sim.body,t0)==CW.motion(s.sim.body,t1) && return
    mesh=CW.transport_mesh(s.sim.body,s.dims,t0,t1);prepared=CW.prepare_projection(mesh)
    baseline=legacy_midpoint_flux(s,mesh,t0,t1)
    revised,mapping=CW.waterlily_face_flux(s,mesh,t0,t1)
    roi=[s.volumes[mesh.cells[mesh.left[e]]]>0 ||
         (mesh.right[e]>0 && s.volumes[mesh.cells[mesh.right[e]]]>0) for e in eachindex(revised)]
    cases=Dict()
    for (label,raw) in [("legacy_midpoint",baseline),("endpoint_paired",revised)]
        q,a=CW.compatible_flux(mesh,raw,t1-t0;prepared)
        cases[label]=merge(a,(roi_relative_l2=norm((q.-raw)[roi])/max(norm(raw[roi]),1e-30),))
    end
    record=merge((t0=t0,t1=t1,cases=cases),mapping);push!(records,record)
    if t1==21. || t1==42.;println(JSON3.write(record));flush(stdout);end
end)
mkpath(dirname(ARGS[1]))
CW.write_json(ARGS[1],Dict("status"=>"SHORTENED_PLATEAU_COUPLING_DIAGNOSTIC",
    "co2_solved"=>false,"inputs"=>c,"elapsed_s"=>elapsed,"records"=>records,
    "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "transport_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","Transport.jl")))),
    "limitations"=>["Identical native airflow for both mappings; no CO2 evolution in this probe.",
        "Real 20-second wall travel, but only one second sealed and one second fully open.",
        "Comparison does not establish operating-cycle accuracy or physical validation."]))
println("Finished mapping diagnostic in ",elapsed," s")
