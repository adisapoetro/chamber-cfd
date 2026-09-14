using ChamberWaterLily

function main(args)
    options=Dict{String,String}()
    length(args)%2==0 || error("Use --input configs/current.toml --output results/new_case [--duration 30]")
    for i in 1:2:length(args)
        args[i] in ["--input","--output","--duration"] || error("Unknown option: $(args[i])")
        haskey(options,args[i]) && error("Repeated option: $(args[i])")
        options[args[i]]=args[i+1]
    end
    haskey(options,"--output") || error("--output is required; choose a fresh folder")
    c=load_config(get(options,"--input",joinpath(@__DIR__,"..","configs","current.toml")))
    duration=parse(Float64,get(options,"--duration",string(c["cycle"]["cycles"]*ChamberWaterLily.period(c))))
    run_airflow(c,abspath(options["--output"]);duration_s=duration)
end

main(ARGS)
