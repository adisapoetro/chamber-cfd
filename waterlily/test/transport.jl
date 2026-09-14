function periodic_mesh(n)
    dims=(n,1,1);cells=collect(CartesianIndices(dims));ids=reshape(collect(1:n),dims)
    CW.TransportMesh(dims,vec(cells),ids,collect(1:n),[mod1(i+1,n) for i in 1:n],
        vec(cells),ones(Int,n),ones(n),ones(n),fill(1/n,n),fill(1/n,n),1/n)
end

function ambient_ended_mesh(n)
    cells=vec(collect(CartesianIndices((n,1,1))));ids=reshape(collect(1:n),(n,1,1))
    # Unit-area interval with ambient reservoirs at x=0 and x=1.
    CW.TransportMesh((n,1,1),cells,ids,vcat(1:n-1,1,n),vcat(2:n,0,0),
        vcat(cells[2:end],cells[[1,end]]),ones(Int,n+1),ones(n+1),ones(n+1),
        fill(1/n,n),fill(1/n,n),1/n)
end

@testset "Stagnant ambient reservoirs diffuse without wind" begin
    errors=Float64[]
    for n in [32,64]
        mesh=ambient_ended_mesh(n);x=((1:n).-0.5)./n
        c=400 .-10sin.(pi.*x);D=.1;t=.1
        steps=ceil(Int,t/(.5/n^2));dt=t/steps;budget=0.
        initial=sum(mesh.old_volume.*c)
        for _ in 1:steps
            c,a=CW.transport_step(mesh,c,zeros(n+1),fill(D,n),zeros(n),dt,400.)
            budget=max(budget,abs(a.budget_error_ppm_m3))
        end
        exact=400 .-10exp(-pi^2*D*t).*sin.(pi.*x)
        push!(errors,sqrt(sum(abs2,c.-exact)/n))
        @test sum(mesh.new_volume.*c)>initial
        @test 390<=minimum(c)<=maximum(c)<=400
        @test budget<1e-10
    end
    @test errors[2]<.3errors[1]
    @test errors[2]<.001
end

@testset "Independent chamber, palm, wind, fan and uptake inputs" begin
    c=load_config(INPUT)
    c["chamber"]["width_m"]=5.;c["chamber"]["height_m"]=5.
    c["palm"]["height_m"]=3.1;c["palm"]["crown_radius_x_m"]=2.
    c["air"]["wind_speed_m_s"]=.8;c["air"]["wind_to_deg"]=45.
    c["fans"]["count"]=5;c["exchange"]["net_co2_umol_s"]=-60.
    validate_config(c);state=build_flow(c);scalar=CW.build_scalar(state)
    @test nominal_volume(c)==150.
    @test length(fan_configuration(c))==5
    @test [f.azimuth for f in fan_configuration(c)]==[0,72,144,216,288]
    @test [f.height for f in fan_configuration(c)]≈[.5,1.5,2.5,3.5,4.5]
    @test sum(scalar.leaf.source)*scalar.leaf.molar_density≈-60.
    advance_flow!(state,.1;after_step=(s,a,z)->CW.advance_scalar!(scalar,s,a,z))
    @test sum(state.volumes.*scalar.concentration)/150≈400-60*.1/(scalar.leaf.molar_density*150) atol=1e-9
end

@testset "Scalar advection, diffusion and bounds" begin
    errors=Float64[];diff_errors=Float64[]
    for n in [32,64]
        mesh=periodic_mesh(n);x=((1:n).-0.5)./n;c=400 .+sin.(2pi.*x)
        dt=.2/n;steps=round(Int,.2/dt);u=.5
        for _ in 1:steps
            c,a=CW.transport_step(mesh,c,fill(u,n),zeros(n),zeros(n),dt,400.)
            @test abs(a.budget_error_ppm_m3)<1e-10
        end
        push!(errors,sqrt(sum(abs2,c.-(400 .+sin.(2pi.*(x.-u*.2))))/n))
        @test 399<=minimum(c)<=maximum(c)<=401
        c=400 .+sin.(2pi.*x);D=.1
        for _ in 1:steps;c,_=CW.transport_step(mesh,c,zeros(n),fill(D,n),zeros(n),dt,400.);end
        push!(diff_errors,sqrt(sum(abs2,c.-(400 .+exp(-4pi^2*D*.2).*sin.(2pi.*x)))/n))
        c,_=CW.transport_step(mesh,fill(400.,n),zeros(n),zeros(n),fill(-1/n,n),2.,400.)
        @test maximum(abs,c.-398)<1e-10
        @test_throws ErrorException CW.transport_step(mesh,fill(1.,n),zeros(n),zeros(n),fill(-10.,n),2.,1.)
    end
    @test errors[2]<.6errors[1]
    @test diff_errors[2]<.6diff_errors[1]
    mesh=periodic_mesh(8)
    @test_throws ErrorException CW.transport_step(mesh,fill(400.,8),collect(1.:8.),zeros(8),zeros(8),.1,400.)
end

@testset "Cut geometry, exact closed capacity and impermeable walls" begin
    c=load_config(INPUT);b=chamber_body(c);dims=(28,20,20)
    mesh=CW.transport_mesh(b,dims,0.,.1);labels,exposed=CW.mesh_components(mesh)
    seed=CartesianIndex(15,11,6);label=labels[mesh.ids[seed]]
    interior=findall(==(label),labels)
    @test !(label in exposed)
    @test sum(mesh.old_volume[interior])≈96. atol=1e-10
    @test minimum(mesh.old_volume)>0
    cold=fill(400.,length(mesh.cells));cold[interior].=390.
    q,a=CW.compatible_flux(mesh,zeros(length(mesh.left)),.1)
    out,z=CW.transport_step(mesh,cold,q,ones(length(cold)),zeros(length(cold)),.1,400.)
    @test maximum(abs,out.-cold)<1e-9
    @test abs(z.budget_error_ppm_m3)<1e-8
    source=zeros(length(cold));source[interior].=-mesh.old_volume[interior]./96
    out,z=CW.transport_step(mesh,fill(400.,length(cold)),q,ones(length(cold)),source,10.,400.)
    @test sum(out[interior].*mesh.new_volume[interior])/96≈400-10/96 atol=1e-10
    @test maximum(abs,out[setdiff(1:length(cold),interior)].-400)<1e-10
end

@testset "Moving walls preserve uniform concentration without resetting cells" begin
    c=load_config(INPUT);b=chamber_body(c);dims=(28,20,20)
    for (t,dt) in [(300.,.5),(300.3,.5),(305.,.5),(319.5,.5),(880.,.5),(899.3,.5),(899.5,.5)]
        mesh=CW.transport_mesh(b,dims,t,t+dt)
        @test sum(mesh.old_volume)≈sum(mesh.new_volume) atol=1e-10
        q,a=CW.compatible_flux(mesh,zeros(length(mesh.left)),dt)
        n=length(mesh.cells)
        result,z=CW.transport_step(mesh,fill(400.,n),q,fill(1.6e-5,n),zeros(n),dt,400.)
        @test maximum(abs,result.-400)<1e-9
        @test a.gcl_max_m3_s<1e-9
        @test abs(z.budget_error_ppm_m3)<1e-7
        # Nonuniform inventories must also conserve and remain in their initial bounds.
        cold=[395+5sin(I[1]) for I in mesh.cells]
        result,z=CW.transport_step(mesh,cold,q,fill(.01,n),zeros(n),dt,400.)
        @test minimum(result)>=min(minimum(cold),400.)-1e-7
        @test maximum(result)<=max(maximum(cold),400.)+1e-7
        @test abs(z.budget_error_ppm_m3)<1e-7
    end
end

@testset "Native velocity and geometry use matching times" begin
    state=build_flow(load_config(INPUT));b=state.sim.body;sim=state.sim
    function blended!(array,t,fluid)
        for I in CartesianIndices(sim.flow.p),k in 1:3
            d,_,v=WaterLily.measure(b,WaterLily.loc(k,I,Float64),t/b.dx)
            mu=WaterLily.μ₀(d,sim.ϵ)
            array[I,k]=mu*fluid[k]+(1-mu)*v[k]
        end
    end
    fluid=SVector(.3,-.2,.1)
    for (t0,t1) in [(300.,300.05),(899.95,900.)]
        mesh=CW.transport_mesh(b,state.dims,t0,t1)
        blended!(sim.flow.u⁰,t0,fluid);blended!(sim.flow.u,t1,fluid)
        q,a=CW.waterlily_face_flux(state,mesh,t0,t1)
        expected=mesh.orientation.*mesh.area.*fluid[mesh.direction]
        @test maximum(abs,q.-expected)<1e-11
        @test a.one_sided_extension_faces>0
        @test 0<a.one_sided_extension_area_fraction<1
    end
    # Independently integrate a linearly changing area times a linear velocity.
    t0=305.1;t1=305.15;mesh=CW.transport_mesh(b,state.dims,t0,t1)
    u0=fluid;u1=fluid.+.07
    blended!(sim.flow.u⁰,t0,u0);blended!(sim.flow.u,t1,u1)
    q,a=CW.waterlily_face_flux(state,mesh,t0,t1)
    boxes0=CW.shell_boxes(b,t0);boxes1=CW.shell_boxes(b,t1)
    expected=map(eachindex(q)) do e
        I=mesh.faces[e];k=mesh.direction[e]
        A0=CW.open_face_area(b,I,k,boxes0);A1=CW.open_face_area(b,I,k,boxes1)
        mesh.orientation[e]*(2A0*u0[k]+A0*u1[k]+A1*u0[k]+2A1*u1[k])/6
    end
    @test a.one_sided_extension_faces==0
    @test maximum(abs,q.-expected)<1e-11
    # Verify the native package keeps the actual previous velocity in u⁰.
    fresh=build_flow(load_config(INPUT));initial=copy(fresh.sim.flow.u)
    advance_flow!(fresh,.02)
    @test fresh.sim.flow.u⁰==initial
end

@testset "Native airflow coupled to leaflet uptake through a full diagnostic cycle" begin
    c=load_config(INPUT)
    c["cycle"]["closed_s"]=1.;c["cycle"]["open_phase_s"]=3.
    c["cycle"]["opening_travel_s"]=1.;c["cycle"]["closing_travel_s"]=1.
    state=build_flow(c);scalar=CW.build_scalar(state)
    @test scalar.leaf.triangles==240
    @test scalar.leaf.total_area_m2≈0.7328095439296075
    @test sum(scalar.leaf.source)*scalar.leaf.molar_density≈c["exchange"]["net_co2_umol_s"]
    for t in 1.:4.
        advance_flow!(state,t;after_step=(s,a,z)->CW.advance_scalar!(scalar,s,a,z))
        @test all(isfinite,scalar.concentration)
        @test maximum(scalar.concentration)<=400+1e-8
        @test abs(scalar.audits[end].cumulative_budget_error_ppm_m3)<1e-7
        if t==1.
            mean=sum(state.volumes.*scalar.concentration)/96
            @test mean≈400+sum(scalar.leaf.source)/96 atol=1e-9
        end
    end
    @test all(a.gcl_max_m3_s<1e-8 for a in scalar.audits)
    @test length(scalar.audits)>4
end
