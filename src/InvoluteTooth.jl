module InvoluteTooth # this submodule provides functions for drawing gear involutes
  using TestItems
  using GLMakie
# GLMakie.closeall()
# GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config

  using Roots
  using LaTeXStrings
  using UnitTypes
  using ..Gears

  rad2deg(rad) = rad*180/π

  """
  calculates the x coordinate of a point along the involute
  Note that the handedness of the involute is determined by th-al: if th-al>0 a right hand involute, th-al<0 is a left hand
  """
  ix(r,al,th) = r*( cos(th) + (th-al)*sin(th) )

  """
  calculates the y coordinate of a point along the involute
  """
  iy(r,al,th) = r*( sin(th) - (th-al)*cos(th) )

  """
  finds thMax to intersect with the line of symmetry gamma 
  Gamma is the angle of the tooth's line of symmetry
  Alpha is the angle of the involutes root on the base circle
  Theta0 is the guess of the angle corresponding to the intersection of the involute and gamma
  """
  function calcThetaIntersect(; gm, al, th0=gm+(gm-al)*2)
    #calculate thMax to intersect with the gamma line
    fxi(thc,al,gm) = cos(thc)^2 + 2*(thc-al)*cos(thc)*sin(thc) + (sin(thc)^2-cos(gm)^2)*(thc-al)^2 - cos(gm)^2
    fyi(thc,al,gm) = sin(thc)^2 - 2*(thc-al)*cos(thc)*sin(thc) + (cos(thc)^2-sin(gm)^2)*(thc-al)^2 - sin(gm)^2
    fxz(thc) = fxi(thc, al, gm)^2 + fyi(thc, al, gm)^2 
    return find_zero(fxz, th0, atol=1e-5) 
  end

  function calcThetaRootIntersection(; rrd, rbd, alb, al0 )
    # fbt(thc) = rrd - rbd*sqrt(1+(thc-btb)^2)
    # thbtr = find_zero(fbt, bto, atol=1e-5) 

    fal(thc) = rrd - rbd*sqrt(1+(thc-alb)^2)
    # return find_zero(fal, alb, atol=1e-5) 
    return find_zero(fal,al0, atol=1e-5) 
  end

  """
  Finds the angle when the involute has risen through toothHeight
  Gamma is the angle of the tooth's line of symmetry
  Alpha is the angle of the involute's root on the base circle
  toothHeight is the desired radial height of the tooth beyond the base circle
  """
  function calcThetaHeight(; r, gm, al, toothHeight)
    htal(th) = r * ( cos(th-gm) + (th-al)*sin(th-gm) - 1 )
    fal(th) = toothHeight - htal(th)
    thal = find_zero(fal, gm + (gm-al)*2)
    return thal
  end
  @testitem "calcThetaHeight" begin
    using UnitTypes
    @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=0.3, toothHeight=20 ), 1.238, atol=1e-3)
    @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=1.27, toothHeight=20 ), 0.550, atol=1e-3)
  end


  function getToothProfilePoints(g::G where G<:Gears.AbstractGear; nPerTooth::Int=100)
    nPerInvolute = convert(Int64,round(nPerTooth/2)) #this may be better defined by physical spacing...

    rod = convert(Float64, convert(Inch, g.outside/2))
    rpd = convert(Float64, convert(Inch, g.pitch/2))
    rrd = convert(Float64, convert(Inch, g.root/2))
    rbd = convert(Float64, convert(Inch, g.base/2))
    t = π/2/24
  
    xs = zeros(g.nTeeth*nPerInvolute*2+1)
    ys = zeros(g.nTeeth*nPerInvolute*2+1)
    for i in 1:round(g.nTeeth)
      gm=2*π/g.nTeeth*(i) #gamma is the angle of the line of the tooth center

      psi = acos(rbd/rpd)
      vpsi = tan(psi)-psi

      # tooth frontside
      alb = gm - asin(t/2/rpd) - vpsi
      alo = Gears.InvoluteTooth.calcThetaHeight(r=rbd, gm=gm, al=alb, toothHeight=rod-rbd) # this is the angle of the involute tip at the outside diameter
      thalr = Gears.InvoluteTooth.calcThetaRootIntersection(rrd=rrd, rbd=rbd, alb=alb, al0=alb)
      ths = LinRange(thalr, alo, nPerInvolute)
      xs[(i-1)*2*nPerInvolute .+ (1:nPerInvolute)] = Gears.InvoluteTooth.ix.(rbd, alb, ths)
      ys[(i-1)*2*nPerInvolute .+ (1:nPerInvolute)] = Gears.InvoluteTooth.iy.(rbd, alb, ths)

      # tooth backside
      btb = gm + asin(t/2/rpd) + vpsi
      bto = Gears.InvoluteTooth.calcThetaHeight(r=rbd, gm=gm, al=btb, toothHeight=rod-rbd) 
      thbtr = Gears.InvoluteTooth.calcThetaRootIntersection(rrd=rrd, rbd=rbd, alb=btb, al0=bto)
      ths = LinRange(bto, thbtr, nPerInvolute)  # reverse direction
      xs[(i-1)*2*nPerInvolute + nPerInvolute .+ (1:nPerInvolute)] = Gears.InvoluteTooth.ix.(rbd, btb, ths)
      ys[(i-1)*2*nPerInvolute + nPerInvolute .+ (1:nPerInvolute)] = Gears.InvoluteTooth.iy.(rbd, btb, ths)

    end
    xs[length(xs)] = xs[1] # close the gear
    ys[length(ys)] = ys[1]

    return (xs,ys)
  end
  @testitem "getToothProfilePoints" begin
    using UnitTypes
    # using GLMakie
    # GLMakie.closeall()
    # GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)
    # fig = Figure()#backgroundcolor="#bbb", size=(1000,1000))
    # axs = Axis(fig[1,1])#, xlabel="X", ylabel="Y", aspect=:data)
    # l = lines!(axs, xs, ys, linewidth=3, label="g")
    # s = scatter!(axs, xs, ys, markersize=8)
    # display(GLMakie.Screen(), fig) # this may not stay open within the test..

    @test isapprox( xs[1], 1.404033, atol=1e-3 ) # values from correct-looking plot
    @test isapprox( xs[2], 1.411619, atol=1e-3 ) 
    @test isapprox( xs[3], 1.419892, atol=1e-3 ) 
    @test isapprox( last(xs), 1.404033, atol=1e-3 ) 
    @test isapprox( ys[1], 0.079214, atol=1e-3 ) # values from correct-looking plot
    @test isapprox( ys[2], 0.081518, atol=1e-3 ) 
    @test isapprox( ys[3], 0.084247, atol=1e-3 ) 
    @test isapprox( last(ys), 0.079214, atol=1e-3 ) 
  end


  function writeToothProfilePoints(g::G where G<:Gears.AbstractGear; fileName::String="", fileExtension::String="txt", nPerTooth::Int=100)
    (xs, ys) = getToothProfilePoints(g, nPerTooth=nPerTooth)

    #now emit the points in the format expected by CAD
    zs = xs *0

    if fileName == ""
      fileName = "gearProfilePoints_$(g.diametral)_$(g.nTeeth).sldcrv"
    end

    #solidworks format: 25mmTAB0mmTAB0mmCRLF
    open(fileName, "w") do f
      for i in 1:length(xs)
        # write(f, @sprintf("%3.3fmm\t%3.3fmm\t%3.3fmm", xs[i],ys[i],zs[i]))
        # print("$(xs[i])mm\t$(ys[i])mm\t$(zs[i])mm\r\n")
        # write(f, "$(xs[i])mm\t$(ys[i])mm\t$(zs[i])mm\r\n")
        write(f, "$(xs[i])in\t$(ys[i])in\t$(zs[i])in\r\n")
      end
    end
  end


  """
  Plots the base arc under the tooth
  """
  function plotBaseSection(; r, al, bt, thExtra=0.5, n=100, axs=nothing, linecolor=:gray, linestyle=:solid)
    if axs == nothing
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    #want to plot only a section of the base circle, that from al to thMax
    ths = LinRange(al-thExtra, bt+thExtra, n) 
    xs = r.*cos.(ths)
    ys = r.*sin.(ths)
    # p = plot(p, xs, ys, linecolor=linecolor, linestyle=linestyle)
    l = lines!(ax, xs, ys, color=linecolor, linestyle=linestyle)

    if axs == nothing
      display(GLMakie.Screen(), fig)
    end
    return l
  end

  """
  Plots the line of symmetry gamma
  """
  function plotGamma(; axs, gm, r, ht, linecolor=:green, linestyle=:dash )
    # p = plot(p, [0,(r+ht)*cos(gm)], [0,(r+ht)*sin(gm)], linecolor=linecolor, linestyle=linestyle)
    l = lines!(axs, [0,(r+ht)*cos(gm)], [0,(r+ht)*sin(gm)], color=linecolor, linestyle=linestyle)
    return l
  end

  """
  Plots the line to the base of the involute
  """
  function plotAlpha(; axs, r, al , linecolor=:gray, linestyle=:dash)
    # p = plot(p, [0,r*cos(al)], [0,r*sin(al)], linecolor=linecolor, linestyle=linestyle)
    l = lines!(axs, [0,r*cos(al)], [0,r*sin(al)], color=linecolor, linestyle=linestyle)
    return l
  end

  """
  Plots an involute of the circle
  """
  function plotInvolute(; r, al, thMax, axs=nothing, linecolor=:blue, linestyle=:solid)
    if axs == nothing
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    ths = LinRange(al, thMax,100) 

    xs = ix.(r, al, ths)
    ys = iy.(r, al, ths)
    # p = plot!(xs,ys, linecolor=linecolor, linestyle=linestyle)
    l = lines!(ax,xs,ys, color=linecolor, linestyle=linestyle)

    if axs == nothing
      display(GLMakie.Screen(), fig)
    end

    return l
  end

  """
  """
  function plotInvoluteConstruction(; r, al, th, axs=nothing, linecolor=:gray, linestyle=:dash)
    if axs == nothing
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    xs = [0, r*cos(th), ix(r,al,th)]
    ys = [0, r*sin(th), iy(r,al,th)]
    # p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle)
    l = lines!(ax, xs,ys, color=linecolor, linestyle=linestyle)

    if axs == nothing
      display(GLMakie.Screen(), fig)
    end

    return l
  end

  """
  """
  function drawAngleArc(; axs, r, aMax, aMin=0, label="label", aLabel=aMax/2, linecolor=:gray, linestyle=:solid, fontsize=10)
    ths = LinRange(aMin, aMax, 100)
    xs = r .* cos.(ths)
    ys = r .* sin.(ths)
    # p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle, markerstyle=:rightarrow)
    s = scatterlines!(axs, xs, ys, color=linecolor, linestyle=linestyle, marker=:rtriangle)

    # annotate!(r*cos(aLabel), r*sin(aLabel), "  "*label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationrotation=rad2deg(aLabel) )
    # annotate!(r*cos(aLabel), r*sin(aLabel), label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationhalign=:left, annotationvalign=:bottom )
    text!(axs, r*cos(aLabel), r*sin(aLabel), text=label, color=linecolor, fontsize=fontsize, align=(:left,:bottom) )
    return s
  end

  """
  """
  function drawToothTop(;r, al, thal, bt, thbt, axs=nothing, linecolor=:orange, linestyle=:dash)
    if axs == nothing
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    # p = plot(p, [ix(r,al,thal), ix(r,bt,thbt)], [iy(r,al,thal), iy(r,bt,thbt)], linecolor=linecolor, linestyle=linestyle)
    l = lines!(ax, [ix(r,al,thal), ix(r,bt,thbt)], [iy(r,al,thal), iy(r,bt,thbt)], color=linecolor, linestyle=linestyle)

    if axs == nothing
      display(GLMakie.Screen(), fig)
    end
    return l
  end

  """
  plots both sides to form a tooth
  r is the root diameter
  gamma is the center angle of the tooth
  delta is the angular width of the tooth at the root diameter
  height is the tooth height, or the difference between the outside and root diameters
  """
  function plotTooth(; r, gm, dl, ht, axs=nothing)
    # plot construction of both involutes
    al = gm - dl/2
    althMax = calcThetaHeight(r=r, gm=gm, al=al, toothHeight=ht)
    bt = gm + dl/2
    btthMax = calcThetaHeight(r=r, gm=gm, al=bt, toothHeight=ht)

    # if axs === nothing
      # p = plot(legend=false, title="Involute tooth with γ=$(round(rad2deg(gm))), δ=$(round(rad2deg(dl)))")
    # end
    if axs == nothing
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end


    lines!(ax, [0,r*1.1],[0,0], color=:black)
    plotBaseSection(axs=ax, r=r, al=al, bt=bt, linecolor=:black, linestyle=:solid)
    plotGamma(axs=ax, r=r, gm=gm, ht=ht)

    # drawAngleArc(axs=ax, r=40, aMin=0, aMax=gm, aLabel=deg2rad(0), label="γ", color=:gray)
    # savefig("on220826_involutes_baseGamma.svg")
    # drawAngleArc(axs=ax, r=70, aMin=al, aMax=bt, aLabel=deg2rad(35), label=L"\delta", color=:gray)
    # savefig("on220826_involutes_delta.svg")

    plotAlpha(axs=ax, r=r, al=al)
    # p = drawAngleArc(p=p, r=50, aMin=0, aMax=al, aLabel=deg2rad(0), label="α", linecolor=:gray)
    plotInvolute(axs=ax, r=r, al=al, thMax=althMax, linecolor=:blue )
    # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(31), linecolor=:blue )
    # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(106), linecolor=:blue )
    plotInvoluteConstruction(axs=ax, r=r, al=al, th=althMax, linecolor=:cyan, linestyle=:solid )
    # p = drawAngleArc(p=p, r=60, aMin=0, aMax=althMax, aLabel=deg2rad(0), label=L"\theta_\alpha", linecolor=:gray)
    # savefig("on220826_involutes_alpha.svg")

    plotAlpha(axs=ax, r=r, al=bt)
    # p = drawAngleArc(p=p, r=30, aMin=0, aMax=bt, aLabel=deg2rad(0), label="β", linecolor=:gray)
    plotInvolute(axs=ax, r=r, al=bt, thMax=btthMax, linecolor=:red )
    plotInvoluteConstruction(axs=ax, r=r, al=bt, th=btthMax, linecolor=:magenta, linestyle=:solid )
    # p = drawAngleArc(p=p, r=70, aMin=0, aMax=btthMax, aLabel=deg2rad(0), label=L"\theta_\beta", linecolor=:gray)
    # savefig("on220826_involutes_beta.svg")

    drawToothTop(axs=ax, r=r, al=al, thal=althMax, bt=bt, thbt=btthMax, linecolor=:orange, linestyle=:solid)
  end



end #InvoluteTooth
