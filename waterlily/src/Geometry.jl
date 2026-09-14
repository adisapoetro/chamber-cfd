"""Two sliding shells. Distances are in grid cells; body velocities are in m/s."""
struct ChamberBody <: WaterLily.AbstractBody
    dx::Float64
    lower::SVector{3,Float64}
    size::SVector{3,Float64}
    floor::Float64
    wall::Float64
    gap::Float64
    timing::NTuple{4,Float64}
end

function chamber_body(c)
    ch=c["chamber"];d=c["domain"];cy=c["cycle"]
    ChamberBody(c["numerics"]["cell_m"],SVector(-d["x_m"]/2,-d["y_m"]/2,0.),
        SVector(ch["width_m"],ch["depth_m"],ch["height_m"]),ch["floor_m"],
        wall_thickness(c),ch["maximum_gap_m"],
        Tuple(Float64.([cy["closed_s"],cy["open_phase_s"],cy["opening_travel_s"],cy["closing_travel_s"]])))
end

@inline function motion(b::ChamberBody,t)
    closed,opened,opening,closing=b.timing
    p=mod(t,closed+opened)
    p<closed && return 0.,0.
    p<closed+opening && return b.gap*(p-closed)/opening,b.gap/opening
    p<closed+opened-closing && return b.gap,0.
    return b.gap*(closed+opened-p)/closing,-b.gap/closing
end

@inline function box_measure(x,center,half)
    y=x-center;q=abs.(y)-half;r=max.(q,0.);mag=norm(r)
    d=mag+min(maximum(q),0.)
    if mag>0
        return d,sign.(y).*r/mag
    end
    axis=argmax(q)
    return d,SVector(ntuple(i->i==axis ? (y[i]>=0 ? 1. : -1.) : 0.,3))
end

@inline function shell_measure(b,x,gap,side)
    W,D,H=b.size;w=b.wall
    outer_lo=side>0 ? gap/2 : -gap/2-W/2-w
    outer_hi=side>0 ? gap/2+W/2+w : -gap/2
    inner_lo=side>0 ? gap/2-2w : -gap/2-W/2
    inner_hi=side>0 ? gap/2+W/2 : -gap/2+2w
    z=b.floor+H/2
    a,na=box_measure(x,SVector((outer_lo+outer_hi)/2,0.,z),SVector((outer_hi-outer_lo)/2,D/2+w,H/2+w))
    d,nd=box_measure(x,SVector((inner_lo+inner_hi)/2,0.,z),SVector((inner_hi-inner_lo)/2,D/2,H/2))
    return a>=-d ? (a,na) : (-d,-nd)
end

@inline function WaterLily.measure(b::ChamberBody,x,t=0.;fastd²=Inf)
    p=b.lower+b.dx*x;gap,speed=motion(b,t*b.dx)
    if gap==0
        center=SVector(0.,0.,b.floor+b.size[3]/2)
        a,na=box_measure(p,center,b.size/2 .+ b.wall)
        d,nd=box_measure(p,center,b.size/2)
        best,n=a>=-d ? (a,na) : (-d,-nd)
        v=zero(x)
    else
        left,nl=shell_measure(b,p,gap,-1)
        right,nr=shell_measure(b,p,gap,1)
        best,n,v=left<=right ? (left,nl,SVector(-speed/2,0.,0.)) : (right,nr,SVector(speed/2,0.,0.))
    end
    # The external ground is z=0, distinct from the raised chamber floor.
    if p[3]<best
        best=p[3];n=SVector(0.,0.,1.);v=zero(x)
    end
    d=best/b.dx
    d^2>fastd² && return d,zero(x),zero(x)
    return d,n,v
end
WaterLily.sdf(b::ChamberBody,x,t=0.;kwargs...)=first(WaterLily.measure(b,x,t;kwargs...))

"""Exact intersection with the fixed nominal chamber region, independent of dx."""
function reference_volumes(c,dims,lower)
    dx=c["numerics"]["cell_m"];ch=c["chamber"]
    lo=SVector(-ch["width_m"]/2,-ch["depth_m"]/2,ch["floor_m"])
    hi=lo+SVector(ch["width_m"],ch["depth_m"],ch["height_m"])
    v=zeros(Float64,dims)
    for I in CartesianIndices(v)
        a=lower+SVector(Tuple(I).-1)*dx;b=a .+ dx
        v[I]=prod(max.(0.,min.(b,hi)-max.(a,lo)))
    end
    return v
end
