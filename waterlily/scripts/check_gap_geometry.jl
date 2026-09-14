# Geometry-only diagnostic of the native diffuse indicator as shells meet.
using ChamberWaterLily, WaterLily, StaticArrays, JSON3, SHA
const CW=ChamberWaterLily
length(ARGS)==1 || error("Give a fresh JSON output path")
ispath(ARGS[1]) && error("Output exists")
cases=Any[]
for dx in [.5,.25]
    c=load_config(joinpath(@__DIR__,"..","configs","current.toml"))
    c["numerics"]["cell_m"]=dx;c["domain"]["z_m"]=10.
    c["waterlily"]["numerical_wall_cells"]=1.5/dx
    b=chamber_body(c);dims=Tuple(round.(Int,[14,10,10]./dx));samples=Any[]
    for gap in [2.,1.,.5,.25,.125,.05,.01,1e-6,0.]
        t=gap==0 ? 0. : CW.period(c)-gap*c["cycle"]["closing_travel_s"]/b.gap
        sharp=sum(CW.cut_volumes(b,dims,t));kernel=0.
        for I in CartesianIndices(dims)
            x=WaterLily.loc(0,I+CartesianIndex(1,1,1),Float64)
            d=WaterLily.sdf(b,x,t/dx)
            kernel+=WaterLily.μ₀(d,c["waterlily"]["kernel_width_cells"])*dx^3
        end
        # Native face kernel at the centre seam, halfway through the thick roof.
        probe=SVector(0.,0.,b.floor+b.size[3]+b.wall/2)
        d=WaterLily.sdf(b,(probe-b.lower)/dx,t/dx)
        push!(samples,Dict("gap_m"=>gap,"sharp_domain_fluid_volume_m3"=>sharp,
            "sampled_domain_kernel_integral_m3"=>kernel,
            "roof_seam_kernel_mobility"=>WaterLily.μ₀(d,c["waterlily"]["kernel_width_cells"])))
    end
    exact=[s["sharp_domain_fluid_volume_m3"] for s in samples]
    maximum(exact)-minimum(exact)<1e-8 || error("Sharp volume changes despite constant displaced volume")
    push!(cases,Dict("cell_m"=>dx,"fixed_numerical_wall_m"=>1.5,"samples"=>samples,
        "near_closed_kernel_jump_m3"=>samples[end-1]["sampled_domain_kernel_integral_m3"]-
            samples[end]["sampled_domain_kernel_integral_m3"]))
end
mkpath(dirname(ARGS[1]))
report=Dict("status"=>"GEOMETRY_DIAGNOSTIC_ONLY","cases"=>cases,
    "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "geometry_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","Geometry.jl")))),
    "interpretation"=>"The sampled kernel integral is a geometry diagnostic, not native fluid mass or physical air capacity.",
    "limitations"=>["No airflow or CO2 solve. Numerical walls and exterior box held fixed.",
        "A persistent kernel jump identifies unresolved interface support; this probe does not qualify a replacement wall treatment."])
CW.write_json(ARGS[1],report)
for c in cases;println(JSON3.write(Dict("cell_m"=>c["cell_m"],"near_closed_kernel_jump_m3"=>c["near_closed_kernel_jump_m3"])));end
