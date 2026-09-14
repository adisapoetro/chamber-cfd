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
            # Ambient on inflow, zero diffusive gradient on outflow/tangential faces.
            g=q[e]<0 ? dt*mesh.area[e]*diffusivity[a]/(mesh.dx/2) : 0.
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
        g=q[e]<0 ? mesh.area[e]*diffusivity[a]/(mesh.dx/2) : 0.
        outward+=q[e]*(q[e]>=0 ? result[a] : ambient)+g*(result[a]-ambient)
    end
    budget=sum(mesh.new_volume.*result)-sum(mesh.old_volume.*cold)-dt*(sum(source)-outward)
    return result,(budget_error_ppm_m3=budget,boundary_outward_ppm_m3_s=outward,
        linear_residual_ppm_m3=residual,source_ppm_m3_s=sum(source))
end

"""Map native BDIM velocity to sharp fluid apertures, then enforce the GCL.

The blended body contribution is removed and divided by native mobility before
integrating over the sharp open area. This deblending is an approximation at
the BDIM interface; both its magnitude and the subsequent correction are saved.
"""
function waterlily_face_flux(state,mesh,t)
    sim=state.sim;b=sim.body;raw=zeros(length(mesh.left))
    for e in eachindex(raw)
        I=mesh.faces[e]+CartesianIndex(1,1,1);k=mesh.direction[e]
        d,_,v=WaterLily.measure(b,WaterLily.loc(k,I,Float64),t/b.dx)
        mu=WaterLily.μ₀(d,sim.ϵ)
        # Time-averaged apertures may be open before a face is covered at step end.
        # Sample mobility/body motion at the midpoint with the averaged velocity.
        mu>1e-8 || error("Open scalar face has no native BDIM support")
        fluid=(sim.flow.u[I,k]-(1-mu)*v[k])/mu
        raw[e]=mesh.orientation[e]*mesh.area[e]*fluid
    end
    return raw
end
