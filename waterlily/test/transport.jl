function periodic_mesh(n)
    dims=(n,1,1);cells=collect(CartesianIndices(dims));ids=reshape(collect(1:n),dims)
    CW.TransportMesh(dims,vec(cells),ids,collect(1:n),[mod1(i+1,n) for i in 1:n],
        vec(cells),ones(Int,n),ones(n),ones(n),fill(1/n,n),fill(1/n,n),1/n)
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
