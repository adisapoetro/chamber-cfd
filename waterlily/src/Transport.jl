"""Backward-Euler conservative upwind advection and face diffusion.

Concentrations are ppm; inventories are ppm m³. This monotone first-order
extension solves small cut cells implicitly without merging or concentration
clipping. Prescribed sinks that exhaust CO₂ are rejected.
"""
function transport_step(mesh,cold,q,diffusivity,source,dt,ambient)
    n=length(mesh.cells);length(cold)==length(source)==n || error("Scalar shape mismatch")
    all(diffusivity.>=0) || error("Negative diffusivity")
    rows=collect(1:n);cols=collect(1:n);vals=copy(mesh.new_volume)
    rhs=mesh.old_volume.*cold.+dt.*source
    gcl=mesh.new_volume.-mesh.old_volume.+dt.*flux_divergence(mesh,q)
    maximum(abs,gcl)<=1e-8 || error("Scalar flux violates the geometric conservation law")
    # Shift by ambient to avoid loss of precision in the dominant 400 ppm baseline.
    rhs .-= mesh.old_volume.*ambient
    for e in eachindex(q)
        a=mesh.left[e];z=mesh.right[e];f=dt*q[e]
        if z>0
            g=dt*mesh.area[e]*(2diffusivity[a]*diffusivity[z]/max(diffusivity[a]+diffusivity[z],eps()))/mesh.dx
            append!(rows,[a,a,z,z]);append!(cols,[a,z,a,z])
            append!(vals,[max(f,0.)+g,min(f,0.)-g,-max(f,0.)-g,-min(f,0.)+g])
        else
            # Ambient on inflow and stagnant far-field faces; zero gradient on outflow.
            # Molecular diffusion remains active even when advection is exactly zero.
            g=q[e]<=0 ? dt*mesh.area[e]*diffusivity[a]/(mesh.dx/2) : 0.
            push!(rows,a);push!(cols,a);push!(vals,max(f,0.)+g)
        end
    end
    matrix=sparse(rows,cols,vals,n,n)
    # Include the measured GCL residual; a uniform field is not imposed by algebra.
    rhs .-= ambient.*gcl
    shifted=matrix\rhs;result=shifted.+ambient
    residual=maximum(abs,matrix*shifted-rhs)
    all(isfinite,result) && minimum(result)>=-1e-9 || error("Prescribed exchange exhausts CO₂ or solve failed")
    residual<=1e-7 || error("Scalar linear solve did not converge")
    outward=0.
    for e in eachindex(q)
        a=mesh.left[e];mesh.right[e]>0 && continue
        g=q[e]<=0 ? mesh.area[e]*diffusivity[a]/(mesh.dx/2) : 0.
        outward+=q[e]*(q[e]>=0 ? result[a] : ambient)+g*(result[a]-ambient)
    end
    budget=sum(mesh.new_volume.*result)-sum(mesh.old_volume.*cold)-dt*(sum(source)-outward)
    return result,(budget_error_ppm_m3=budget,boundary_outward_ppm_m3_s=outward,
        linear_residual_ppm_m3=residual,source_ppm_m3_s=sum(source))
end

"""Map native BDIM velocity to time-integrated sharp fluid apertures.

During motion, each endpoint velocity is paired with its own geometry and body
velocity. Linear interpolation of supported endpoint fluid velocities is
integrated against the moving area with two-point Gauss quadrature on each
linear geometry interval. A newly exposed/covered face uses the one supported
endpoint; that first-order extension is counted, never replaced by ambient flow.
Static steps use the end velocity, as in backward-Euler scalar transport.

Deblending omits the BDIM first-moment term, so it remains an approximation.
The subsequent conservative projection and its failure screen are still needed.
"""
function waterlily_face_flux(state,mesh,t0,t1)
    sim=state.sim;b=sim.body;raw=zeros(length(mesh.left))
    moving=motion(b,t0)!=motion(b,t1)
    samples=Tuple{Float64,Float64,Vector{RectSolid}}[]
    if moving
        knots=geometry_times(b,state.dims,t0,t1)
        for j in 1:length(knots)-1
            mid=(knots[j]+knots[j+1])/2;half=(knots[j+1]-knots[j])/2
            for g in [-inv(sqrt(3.)),inv(sqrt(3.))]
                t=mid+half*g
                push!(samples,(half/(t1-t0),(t-t0)/(t1-t0),shell_boxes(b,t)))
            end
        end
    end
    extension_faces=0;extension_area=0.
    for e in eachindex(raw)
        I=mesh.faces[e]+CartesianIndex(1,1,1);k=mesh.direction[e]
        x=WaterLily.loc(k,I,Float64)
        d1,_,v1=WaterLily.measure(b,x,t1/b.dx);m1=WaterLily.μ₀(d1,sim.ϵ)
        if !moving
            m1>1e-8 || error("Open scalar face has no native BDIM support")
            raw[e]=mesh.orientation[e]*mesh.area[e]*(sim.flow.u[I,k]-(1-m1)*v1[k])/m1
            continue
        end
        d0,_,v0=WaterLily.measure(b,x,t0/b.dx);m0=WaterLily.μ₀(d0,sim.ϵ)
        m0>1e-8 || m1>1e-8 || error("Open moving scalar face lacks support at both endpoints")
        u0=m0>1e-8 ? (sim.flow.u⁰[I,k]-(1-m0)*v0[k])/m0 : NaN
        u1=m1>1e-8 ? (sim.flow.u[I,k]-(1-m1)*v1[k])/m1 : NaN
        if !isfinite(u0) || !isfinite(u1)
            extension_faces+=1;extension_area+=mesh.area[e]
            isfinite(u0) ? (u1=u0) : (u0=u1)
        end
        raw[e]=mesh.orientation[e]*sum(w*open_face_area(b,mesh.faces[e],k,boxes)*((1-alpha)*u0+alpha*u1)
                                       for (w,alpha,boxes) in samples)
    end
    return raw,(one_sided_extension_faces=extension_faces,
                one_sided_extension_area_fraction=extension_area/sum(mesh.area))
end
