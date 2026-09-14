function clip_polygon(poly,k,bound,above)
    isempty(poly) && return poly
    out=SVector{3,Float64}[];a=poly[end];ina=above ? a[k]>=bound : a[k]<=bound
    for z in poly
        inz=above ? z[k]>=bound : z[k]<=bound
        if ina!=inz
            push!(out,a+(z-a)*((bound-a[k])/(z[k]-a[k])))
        end
        inz && push!(out,z)
        a=z;ina=inz
    end
    return out
end

function polygon_area(poly)
    length(poly)<3 && return 0.
    return sum(norm(cross(poly[j]-poly[1],poly[j+1]-poly[1]))/2 for j in 2:length(poly)-1)
end

"""Distribute the prescribed net exchange by leaflet area intersecting each cell.

Only source-tagged triangles contribute. Trunk and bare frond illustrations
have no source. The geometry is representative, not measured leaf area.
"""
function leaflet_source(state)
    c=state.c;b=state.sim.body;dx=b.dx
    raw=JSON3.read(read(joinpath(@__DIR__,"..","assets","reference_palm.json"),String))
    env=raw.reference_envelope
    scale=SVector(c["palm"]["crown_radius_x_m"]/env.crown_radii_m[1],
        c["palm"]["crown_radius_y_m"]/env.crown_radii_m[2],c["palm"]["height_m"]/env.height_m)
    area=zeros(state.dims);total=0.;triangles=0
    for (i,triangle) in enumerate(raw.triangles)
        raw.source_mask[i] || continue
        poly=[(SVector{3,Float64}(p)-SVector(0.,0.,env.floor_m)).*scale+SVector(0.,0.,b.floor) for p in triangle]
        total+=polygon_area(poly);triangles+=1
        lo=reduce((a,z)->min.(a,z),poly);hi=reduce((a,z)->max.(a,z),poly)
        firstcell=max.(1,floor.(Int,(lo-b.lower)/dx).+1)
        lastcell=min.(SVector(state.dims),floor.(Int,(hi-b.lower)/dx).+1)
        for I in CartesianIndices(Tuple(firstcell[k]:lastcell[k] for k in 1:3))
            cut=poly;a=b.lower+SVector(Tuple(I).-1)*dx
            for k in 1:3
                cut=clip_polygon(cut,k,a[k],true);cut=clip_polygon(cut,k,a[k]+dx,false)
            end
            area[I]+=polygon_area(cut)
        end
    end
    isapprox(sum(area),total;rtol=1e-10) || error("Leaflet area falls outside the domain")
    total>0 || error("Empty leaflet source")
    # µmol/s divided by mol/m³ equals ppm m³/s.
    molar_density=c["air"]["pressure_pa"]/(8.31446261815324*c["air"]["temperature_k"])
    source=area.*(c["exchange"]["net_co2_umol_s"]/(molar_density*total))
    return (;source,area,total_area_m2=total,triangles,molar_density)
end
