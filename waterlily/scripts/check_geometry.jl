using ChamberWaterLily, WaterLily, StaticArrays, JSON3

"""Probe the closed native face-mobility component; no scalar is evolved."""
function geometry_probe(dx)
    c=load_config(joinpath(@__DIR__,"..","configs","current.toml"))
    c["numerics"]["cell_m"]=dx
    state=build_flow(c);sim=state.sim;b=sim.body
    alpha=zeros(Float64,size(sim.flow.p))
    for I in WaterLily.inside(alpha)
        d=WaterLily.sdf(b,WaterLily.loc(0,I,Float64),0.)
        alpha[I]=WaterLily.μ₀(d,sim.ϵ)
    end
    center=SVector(0.,0.,c["chamber"]["floor_m"]+c["chamber"]["height_m"]/2)
    seed=CartesianIndex(Tuple(floor.(Int,(center-b.lower)/dx).+2))
    queue=[seed];seen=falses(size(alpha));seen[seed]=true;next=1
    while next<=length(queue)
        I=queue[next];next+=1
        for k in 1:3,sign in [-1,1]
            J=I+sign*WaterLily.δ(k,I)
            all(2<=J[d]<=size(alpha,d)-1 for d in 1:3) || continue
            face=sign>0 ? J : I
            if !seen[J] && alpha[J]>0 && sim.flow.μ₀[face,k]>0
                seen[J]=true;push!(queue,J)
            end
        end
    end
    reaches_boundary=any(any((I[d]==2 && sim.flow.μ₀[I,d]>0) ||
        (I[d]==size(alpha,d)-1 && sim.flow.μ₀[I+WaterLily.δ(d,I),d]>0) for d in 1:3) for I in queue)
    capacity=sum(alpha[I] for I in queue)*dx^3
    return Dict("cell_m"=>dx,"numerical_wall_m"=>ChamberWaterLily.wall_thickness(c),
        "reference_volume_m3"=>sum(state.volumes),"sampled_kernel_capacity_m3"=>capacity,
        "relative_kernel_capacity_difference"=>capacity/nominal_volume(c)-1,
        "closed_component_reaches_domain_boundary"=>reaches_boundary,
        "kernel_capacity_interpretation"=>"Postprocessed candidate scalar capacity; not a native WaterLily CO2 state or validated air volume.")
end

length(ARGS)==1 || error("Give a new output JSON path")
ispath(ARGS[1]) && error("Output exists")
cases=[geometry_probe(dx) for dx in [.5,.4,.25]]
report=Dict("scientific_status"=>"DIAGNOSTIC_NOT_VALIDATED","co2_transport_enabled"=>false,"cases"=>cases)
mkpath(dirname(abspath(ARGS[1])))
ChamberWaterLily.write_json(ARGS[1],report)
JSON3.pretty(stdout,report);println()
