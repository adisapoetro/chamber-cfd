smagorinsky(I; S, Cs, Δ) = (Cs*Δ)^2*sqrt(2sum(abs2,@view S[I,:,:]))

"""Prescribed momentum sources passed through WaterLily's native udf interface."""
function forcing_fields(c,sim)
    b=sim.body;dx=b.dx;shape=size(sim.flow.u)
    fans=fan_configuration(c);fan_field=zeros(Float64,shape)
    drag=zeros(Float64,size(sim.flow.p));S=zeros(Float64,size(sim.flow.p)...,3,3)
    rho=c["air"]["pressure_pa"]/(287.05*c["air"]["temperature_k"])
    ch=c["chamber"];palm=c["palm"];num=c["numerics"]
    Cs=Float64(num["smagorinsky"])
    for fan in fans, k in 1:3
        fan.direction[k]==0 && continue
        weights=zeros(Float64,size(sim.flow.p))
        for I in WaterLily.inside(sim.flow.p)
            x=b.lower+dx*WaterLily.loc(k,I,Float64)
            if abs(x[1])<ch["width_m"]/2 && abs(x[2])<ch["depth_m"]/2 && ch["floor_m"]<x[3]<ch["floor_m"]+ch["height_m"]
                weights[I]=exp(-sum(abs2,x-fan.position)/(2c["fans"]["kernel_sigma_m"]^2))
            end
        end
        total=sum(weights)*dx^3
        total>0 || error("No grid support for fan force")
        @views fan_field[:,:,:,k] .+= weights .* (dx*fan.thrust*fan.direction[k]/(rho*total))
    end
    for I in WaterLily.inside(sim.flow.p)
        x,y,z=b.lower+dx*WaterLily.loc(0,I,Float64);z-=ch["floor_m"]
        h=palm["height_m"];rx=palm["crown_radius_x_m"];ry=palm["crown_radius_y_m"]
        crown=abs(x)<rx && abs(y)<ry && 0.15h<z<h ? exp(-((x/(0.8rx))^4+(y/(0.8ry))^4+((z-0.62h)/(0.32h))^4)) : 0.
        trunk=0<z<0.45h ? exp(-((x/palm["trunk_radius_m"])^2+(y/palm["trunk_radius_m"])^2)) : 0.
        drag[I]=(num["canopy_drag_per_m"]*crown+num["trunk_drag_per_m"]*trunk)*dx
    end
    drag_force=zeros(Float64,shape)
    active=Ref(true);max_rate=Ref(0.);max_eddy=Ref(0.)
    function force!(flow,u,t)
        if Cs>0
            WaterLily.sgs!(flow,u,t;νₜ=smagorinsky,S,Cs,Δ=1.)
        end
        fill!(drag_force,0.);max_rate[]=0.;max_eddy[]=0.
        for I in WaterLily.inside(flow.p)
            velocity=SVector(ntuple(k->(u[I,k]+u[I+WaterLily.δ(k,I),k])/2,3))
            rate=drag[I]*norm(velocity)
            max_rate[]=max(max_rate[],rate)
            max_eddy[]=max(max_eddy[],smagorinsky(I;S,Cs,Δ=1.))
            for k in 1:3
                drag_force[I,k]=-rate*velocity[k]
            end
        end
        for I in WaterLily.inside(flow.p), k in 1:3
            flow.f[I,k]+=(drag_force[I,k]+drag_force[I-WaterLily.δ(k,I),k])/2
            active[] && (flow.f[I,k]+=fan_field[I,k])
        end
    end
    return (;force!,fan_field,drag,S,active,max_rate,max_eddy)
end
