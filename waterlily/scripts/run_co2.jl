using ChamberWaterLily

function main(args)
    input=joinpath(@__DIR__,"..","configs","current.toml");output=nothing;duration=nothing
    i=1
    while i<=length(args)
        key=args[i];i+=1;i<=length(args) || error("Missing value for $key")
        value=args[i];i+=1
        key=="--input" ? (input=value) : key=="--output" ? (output=value) :
            key=="--duration" ? (duration=parse(Float64,value)) : error("Unknown option: $key")
    end
    isnothing(output) && error("Provide --output with a fresh folder")
    c=load_config(input)
    isnothing(duration) ? run_coupled(c,output) : run_coupled(c,output;duration_s=duration)
end
main(ARGS)
