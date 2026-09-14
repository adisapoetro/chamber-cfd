"""Exact axis-aligned fluid intersections for the same thick shells used by BDIM.

The scalar has its own finite-volume geometry; it is not WaterLily's kernel
indicator. Interior dimensions remain exact at every grid spacing.
"""
struct RectSolid
    lo::SVector{3,Float64}
    hi::SVector{3,Float64}
    sign::Float64
end

function shell_boxes(b::ChamberBody,t)
    gap,_=motion(b,t);W,D,H=b.size;w=b.wall;f=b.floor
    boxes=RectSolid[]
    for side in [-1,1]
        lo=side<0 ? -gap/2-W/2-w : gap/2
        hi=side<0 ? -gap/2 : gap/2+W/2+w
        ilo=side<0 ? -gap/2-W/2 : gap/2
        ihi=side<0 ? -gap/2 : gap/2+W/2
        push!(boxes,RectSolid(SVector(lo,-D/2-w,f-w),SVector(hi,D/2+w,f+H+w),1.))
        push!(boxes,RectSolid(SVector(ilo,-D/2,f),SVector(ihi,D/2,f+H),-1.))
    end
    return boxes
end

@inline box_overlap(lo,hi,box)=prod(max.(0.,min.(hi,box.hi)-max.(lo,box.lo)))

function cut_volumes(b,dims,t)
    boxes=shell_boxes(b,t);vol=zeros(Float64,dims);dx=b.dx
    for I in CartesianIndices(vol)
        lo=b.lower+SVector(Tuple(I).-1)*dx;hi=lo .+ dx
        lo=setindex(lo,max(lo[3],0.),3)
        v=prod(max.(0.,hi-lo))-sum(x.sign*box_overlap(lo,hi,x) for x in boxes)
        -1e-12*dx^3<=v<=dx^3*(1+1e-12) || error("Invalid geometric cell volume")
        vol[I]=clamp(v,0.,dx^3) # roundoff of exact rectangular intersections only
        vol[I]<1e-13*dx^3 && (vol[I]=0.)
    end
    return vol
end

function open_face_area(b,I,k,boxes)
    dx=b.dx;lo=b.lower+SVector(Tuple(I).-1)*dx;hi=lo .+ dx
    k==3 && lo[3]<=0 && return 0. # impermeable external ground
    lo=setindex(lo,max(lo[3],0.),3)
    area=dx^2
    # Minimum of both one-sided limits closes a face lying exactly on a wall.
    for side in [-1.,1.]
        solid=0.;x=lo[k]+side*1e-10*dx
        for box in boxes
            if box.lo[k]<=x<box.hi[k]
                solid+=box.sign*prod(max(0.,min(hi[d],box.hi[d])-max(lo[d],box.lo[d])) for d in 1:3 if d!=k)
            end
        end
        area=min(area,dx^2-solid)
    end
    -1e-9*dx^2<=area<=dx^2*(1+1e-9) || error("Invalid face aperture")
    return area<1e-12*dx^2 ? 0. : min(area,dx^2)
end

"""Split a linear shell displacement wherever an x edge crosses a grid face."""
function geometry_times(b,dims,t0,t1)
    t1>t0 || error("Positive geometry interval required")
    start=shell_boxes(b,t0);finish=shell_boxes(b,t1)
    knots=Float64[t0,t1]
    for (a,z) in zip(start,finish),edge in [:lo,:hi]
        x0=getproperty(a,edge)[1];x1=getproperty(z,edge)[1]
        x1==x0 && continue
        for j in 0:dims[1]
            face=b.lower[1]+j*b.dx;s=(face-x0)/(x1-x0)
            1e-10<s<1-1e-10 && push!(knots,t0+s*(t1-t0))
        end
    end
    return sort!(unique!(knots))
end

struct TransportMesh
    shape::NTuple{3,Int}
    cells::Vector{CartesianIndex{3}}
    ids::Array{Int,3}
    left::Vector{Int}
    right::Vector{Int} # zero is the ambient domain boundary
    faces::Vector{CartesianIndex{3}}
    direction::Vector{Int}
    orientation::Vector{Float64}
    area::Vector{Float64}
    old_volume::Vector{Float64}
    new_volume::Vector{Float64}
    dx::Float64
end

function transport_mesh(b,dims,t0,t1)
    old=cut_volumes(b,dims,t0);new=cut_volumes(b,dims,t1)
    cells=findall((old.+new).>0);ids=zeros(Int,dims)
    for (i,I) in enumerate(cells);ids[I]=i;end
    knots=geometry_times(b,dims,t0,t1)
    samples=[((knots[j+1]-knots[j])/(t1-t0),shell_boxes(b,(knots[j+1]+knots[j])/2)) for j in 1:length(knots)-1]
    left=Int[];right=Int[];faces=CartesianIndex{3}[];direction=Int[];orientation=Float64[];areas=Float64[]
    function add(a,z,I,k,sign)
        area=sum(weight*open_face_area(b,I,k,boxes) for (weight,boxes) in samples)
        area<=0 && return
        push!(left,a);push!(right,z);push!(faces,I);push!(direction,k);push!(orientation,sign);push!(areas,area)
    end
    for (a,I) in enumerate(cells),k in 1:3
        if I[k]==1;add(a,0,I,k,-1.);end
        J=I+CartesianIndex(ntuple(d->d==k ? 1 : 0,3))
        if I[k]==dims[k]
            add(a,0,J,k,1.)
        elseif ids[J]>0
            add(a,ids[J],J,k,1.)
        end
    end
    return TransportMesh(dims,cells,ids,left,right,faces,direction,orientation,areas,old[cells],new[cells],b.dx)
end

"""Conservative face incidence: positive flux leaves the left cell."""
function flux_divergence(mesh,q)
    d=zeros(length(mesh.cells))
    for e in eachindex(q)
        a=mesh.left[e];z=mesh.right[e];d[a]+=q[e]
        z>0 && (d[z]-=q[e])
    end
    return d
end

function mesh_components(mesh)
    n=length(mesh.cells);parent=collect(1:n)
    function root(i)
        while parent[i]!=i;parent[i]=parent[parent[i]];i=parent[i];end
        return i
    end
    for e in eachindex(mesh.left)
        mesh.right[e]>0 || continue
        a=root(mesh.left[e]);z=root(mesh.right[e]);parent[a]=z
    end
    labels=[root(i) for i in 1:n]
    exposed=Set(labels[mesh.left[e]] for e in eachindex(mesh.left) if mesh.right[e]==0)
    return labels,exposed
end

"""Least-squares face-flux correction enforcing the discrete geometric law.

It changes scalar fluxes only; momentum and pressure remain native WaterLily.
Neumann components are checked before fixing a gauge; mass is never repaired.
"""
function prepare_projection(mesh)
    n=length(mesh.cells);rows=Int[];cols=Int[];vals=Float64[]
    mobility=mesh.area./mesh.dx
    for e in eachindex(mesh.left)
        a=mesh.left[e];z=mesh.right[e];w=mobility[e]
        push!(rows,a);push!(cols,a);push!(vals,w)
        if z>0
            append!(rows,[z,a,z]);append!(cols,[z,z,a]);append!(vals,[w,-w,-w])
        end
    end
    labels,exposed=mesh_components(mesh);sealed=Vector{Int}[]
    for label in unique(labels)
        label in exposed && continue
        members=findall(==(label),labels);push!(sealed,members)
        # A diagonal gauge is equivalent to p=0 at one cell when compatibility holds.
        push!(rows,members[1]);push!(cols,members[1]);push!(vals,1.)
    end
    L=sparse(rows,cols,vals,n,n)
    return (;factor=cholesky(Symmetric(L)),mobility,sealed)
end

function compatible_flux(mesh,raw,dt;prepared=prepare_projection(mesh))
    target=(mesh.old_volume.-mesh.new_volume)./dt
    rhs=target.-flux_divergence(mesh,raw);compat=0.
    for members in prepared.sealed
        r=sum(rhs[members]);compat=max(compat,abs(r))
        abs(r)<=1e-9*max(1.,sum(abs,raw)) || error("Incompatible sealed component: $r m³/s")
    end
    potential=prepared.factor\rhs
    mobility=prepared.mobility
    corrected=copy(raw)
    for e in eachindex(raw)
        a=mesh.left[e];z=mesh.right[e]
        corrected[e]+=mobility[e]*(potential[a]-(z>0 ? potential[z] : 0.))
    end
    residual=maximum(abs,flux_divergence(mesh,corrected).-target)
    residual<=1e-8 || error("Geometric conservation residual too large: $residual")
    delta=corrected.-raw
    audit=(gcl_max_m3_s=residual,component_compatibility_m3_s=compat,
        correction_relative_l2=norm(delta)/max(norm(raw),1e-30),
        correction_max_m_s=maximum(abs.(delta)./mesh.area),
        raw_max_m_s=maximum(abs.(raw)./mesh.area))
    return corrected,audit
end
