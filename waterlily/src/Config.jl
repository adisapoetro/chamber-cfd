"""Read SI inputs. Missing optional fan lists request automatic placement."""
function load_config(path)
    c = TOML.parsefile(path)
    validate_config(c)
    return c
end

function validate_config(c)
    c["schema_version"] == 1 || error("Unsupported schema_version")
    expected = TOML.parsefile(joinpath(@__DIR__, "..", "configs", "current.toml"))
    Set(keys(c)) == Set(keys(expected)) || error("Missing or unknown configuration group")
    for (group, values) in expected
        values isa AbstractDict || continue
        allowed = union(Set(keys(values)), group == "fans" ? Set(["heights_m", "azimuths_deg"]) : Set{String}())
        all(k in allowed for k in keys(c[group])) || error("Unknown input in $group")
        all(haskey(c[group], k) for k in keys(values)) || error("Missing input in $group")
        for (k, v) in c[group]
            v isa Number && !isfinite(v) && error("Non-finite $group.$k")
            v isa AbstractVector && !all(isfinite, v) && error("Non-finite $group.$k")
        end
    end
    for (g, names) in [("chamber", ["width_m", "depth_m", "height_m", "maximum_gap_m", "physical_wall_thickness_m"]),
                       ("domain", ["x_m", "y_m", "z_m"]),
                       ("air", ["temperature_k", "pressure_pa"]),
                       ("cycle", ["closed_s", "open_phase_s", "opening_travel_s", "closing_travel_s"]),
                       ("palm", ["height_m", "crown_radius_x_m", "crown_radius_y_m", "trunk_radius_m"]),
                       ("fans", ["diameter_m", "kernel_sigma_m"]),
                       ("numerics", ["cell_m", "dt_s", "save_every_s", "turbulent_schmidt"]),
                       ("waterlily", ["kernel_width_cells", "numerical_wall_cells", "pressure_tolerance"]) ]
        all(c[g][k] > 0 for k in names) || error("Expected positive values in $g")
    end
    for (g, names) in [("chamber", ["floor_m"]), ("air", ["wind_speed_m_s", "ambient_co2_ppm", "initial_co2_ppm"]),
                       ("fans", ["free_air_cfm_per_fan", "momentum_factor"]),
                       ("numerics", ["viscosity_m2_s", "molecular_diffusivity_m2_s", "smagorinsky", "canopy_drag_per_m", "trunk_drag_per_m"])]
        all(c[g][k] >= 0 for k in names) || error("Expected nonnegative values in $g")
    end
    n=c["fans"]["count"]; n isa Integer && n>=0 || error("fans.count must be a nonnegative integer")
    cyc=c["cycle"]; cyc["cycles"] isa Integer && cyc["cycles"]>0 || error("cycle.cycles must be a positive integer")
    cyc["opening_travel_s"]+cyc["closing_travel_s"] <= cyc["open_phase_s"] || error("Travel exceeds open phase")
    c["fans"]["schedule"] == "closed_only" || error("Only closed_only fan operation is implemented")
    c["exchange"]["model"] == "constant_leaflet_net_exchange" || error("Unknown exchange input model")
    c["waterlily"]["numerical_wall_cells"] >= 2c["waterlily"]["kernel_width_cells"]+1 || error("BDIM walls need an exact-zero core across grid phases")
    c["waterlily"]["pressure_max_iterations"] isa Integer && c["waterlily"]["pressure_max_iterations"]>0 || error("Invalid pressure_max_iterations")
    c["display"]["co2_min_ppm"] < c["display"]["co2_max_ppm"] || error("Invalid CO2 colour range")
    c["display"]["fps"] isa Integer && c["display"]["fps"]>0 || error("Invalid display.fps")
    ch=c["chamber"]; p=c["palm"]; w=wall_thickness(c)
    p["height_m"]<ch["height_m"] && 2p["crown_radius_x_m"]<ch["width_m"] && 2p["crown_radius_y_m"]<ch["depth_m"] || error("Palm envelope does not fit the chamber")
    d=c["domain"]
    d["x_m"]>ch["width_m"]+ch["maximum_gap_m"]+2w && d["y_m"]>ch["depth_m"]+2w && d["z_m"]>ch["floor_m"]+ch["height_m"]+w || error("Exterior domain does not contain the numerically thickened moving walls")
    fan_configuration(c)
    return c
end

nominal_volume(c) = prod(c["chamber"][k] for k in ["width_m","depth_m","height_m"])
wall_thickness(c) = max(c["chamber"]["physical_wall_thickness_m"], c["waterlily"]["numerical_wall_cells"]*c["numerics"]["cell_m"])
period(c) = c["cycle"]["closed_s"] + c["cycle"]["open_phase_s"]

function cycle_state(c, t)
    cy=c["cycle"]; p=mod(t,period(c)); g=c["chamber"]["maximum_gap_m"]
    if p<cy["closed_s"]
        return (gap=0.0, gap_speed=0.0, fans_on=true, phase="sealed")
    elseif p<cy["closed_s"]+cy["opening_travel_s"]
        return (gap=g*(p-cy["closed_s"])/cy["opening_travel_s"], gap_speed=g/cy["opening_travel_s"], fans_on=false, phase="opening")
    elseif p<period(c)-cy["closing_travel_s"]
        return (gap=g, gap_speed=0.0, fans_on=false, phase="open")
    end
    return (gap=g*(period(c)-p)/cy["closing_travel_s"], gap_speed=-g/cy["closing_travel_s"], fans_on=false, phase="closing")
end

function next_event(c,t)
    start=floor(t/period(c))*period(c)
    cy=c["cycle"]
    events=start .+ [cy["closed_s"],cy["closed_s"]+cy["opening_travel_s"],period(c)-cy["closing_travel_s"],period(c)]
    return first(e for e in events if e>t+1e-9)
end

function fan_configuration(c)
    f=c["fans"];ch=c["chamber"];n=f["count"]
    heights=get(f,"heights_m",[(i-0.5)*ch["height_m"]/n for i in 1:n])
    angles=get(f,"azimuths_deg",[f["first_azimuth_deg"]+360(i-1)/n for i in 1:n])
    length(heights)==length(angles)==n || error("Fan lists must match fans.count")
    abs(f["pole_x_m"])+f["diameter_m"]/2<ch["width_m"]/2 && abs(f["pole_y_m"])+f["diameter_m"]/2<ch["depth_m"]/2 || error("Fan support lies outside the chamber")
    all(f["diameter_m"]/2<h<ch["height_m"]-f["diameter_m"]/2 for h in heights) || error("Fan disc does not fit vertically")
    rho=c["air"]["pressure_pa"]/(287.05*c["air"]["temperature_k"])
    q=f["free_air_cfm_per_fan"]*0.0004719474432
    thrust=rho*q^2/(pi*f["diameter_m"]^2/4)*f["momentum_factor"]
    return [(position=SVector(f["pole_x_m"],f["pole_y_m"],ch["floor_m"]+Float64(heights[i])),
             height=Float64(heights[i]),azimuth=Float64(angles[i]),
             direction=SVector(cosd(angles[i]),sind(angles[i]),0.0),thrust=thrust) for i in 1:n]
end
