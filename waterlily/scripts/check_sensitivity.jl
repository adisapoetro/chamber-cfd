# Short matched-domain checks; not convergence evidence for the full cycle.
using ChamberWaterLily, WaterLily, JSON3
const CW=ChamberWaterLily

function trial(label;dx=.5,dt=.5,duration=10.,wind=0.,fans=3)
    c=load_config(joinpath(@__DIR__,"..","configs","current.toml"))
    c["numerics"]["cell_m"]=dx;c["numerics"]["dt_s"]=dt
    c["domain"]["z_m"]=10. # Same actual box for both grids.
    c["waterlily"]["numerical_wall_cells"]=1.5/dx # Hold solid geometry fixed.
    c["air"]["wind_to_deg"]=wind;c["fans"]["count"]=fans
    state=build_flow(c);scalar=CW.build_scalar(state)
    elapsed=@elapsed advance_flow!(state,duration;after_step=(s,a,z)->CW.advance_scalar!(scalar,s,a,z))
    mean=sum(state.volumes.*scalar.concentration)/96
    sd=sqrt(sum(state.volumes.*(scalar.concentration.-mean).^2)/96)
    return Dict("case"=>label,"cell_m"=>dx,"dt_max_s"=>dt,"duration_s"=>duration,
        "wind_to_deg"=>wind,"fans"=>fans,"actual_domain_m"=>collect(state.dims).*dx,
        "numerical_wall_m"=>CW.wall_thickness(c),"co2_mean_ppm"=>mean,"co2_spatial_sd_ppm"=>sd,
        "sealed_mean_error_ppm"=>mean-(400+sum(scalar.leaf.source)*duration/96),
        "budget_error_max_ppm_m3"=>maximum(abs(a.cumulative_budget_error_ppm_m3) for a in scalar.audits),
        "roi_correction_max_relative_l2"=>maximum(a.roi_correction_relative_l2 for a in scalar.audits),
        "roi_correction_final_relative_l2"=>scalar.audits[end].roi_correction_relative_l2,
        "steps"=>length(scalar.audits),"elapsed_s"=>elapsed)
end

length(ARGS)==1 || error("Give a fresh JSON output path")
ispath(ARGS[1]) && error("Output exists")
cases=Any[]
for (name,kwargs) in [("reference",(;)),("smaller_time_step",(;dt=.025)),
        ("refined_grid_fixed_wall",(;dx=.25)),("crosswind_no_fans",(;wind=90.,fans=0,duration=1.))]
    r=trial(name;kwargs...);push!(cases,r);println(JSON3.write(r));flush(stdout)
    CW.write_json(ARGS[1],Dict("scientific_status"=>"SHORT_DIAGNOSTIC_NOT_FULL_CYCLE_CONVERGENCE","cases"=>cases))
end
