module Dev
  using UnitTypes
  using Gears
  using GLMakie
  using Roots

  GLMakie.closeall()
  GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config


  function testToothPlotting()
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    # (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)

    fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
    axs = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    # lines!(axs, [0,1],[0,1])

    r = 1
    al = 1
    bt = 2
    gm = 1.5
    ht = 0.5
    thExtra = 0.5

    # Gears.InvoluteTooth.plotTooth(axs=axs, bd=BaseDiameter(Centimeter(2)), gm=Radian(gm), dl=Radian(bt-al), ht=Millimeter(ht) )
    Gears.InvoluteTooth.plotTooth(axs=axs, bd=g.base, gm=Radian(gm), dl=Radian(bt-al), ht=Millimeter(ht) )
    DataInspector(fig, transparency=true, backgroundcolor=RGBAf(1,1,1,0.5)) # https://docs.makie.org/stable/explanations/inspector/index.html
    display(GLMakie.Screen(), fig) # note the window only lasts as long as the julia session
  end
  # testToothPlotting()

  function dev1129()
    g = GearANSI( PitchDiameter(Inch(1.2500)), 30, Degree(20) ) # sdpsi_s10c9z-024h030

    # @show gearUnit = typeof(g.pitch.measure) # for consistency convert every length into this unit
    # @show rpd = g.pitch.measure/2 
    # @show t = gearUnit(π/(2*g.nTeeth))  # calculating the tooth width t at the pitch diameter = pi/(2*Nteeth) cf Dooner61, units are the same as pitch diameter
    # @show t/2/rpd
    # @show Radian(asin((t/2/rpd).value)) 


    (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g)
    fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
    axs = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())

    scatterlines!(axs, toBaseFloat.(xs), toBaseFloat.(ys), linewidth=1)

    DataInspector(fig, transparency=true, backgroundcolor=RGBAf(1,1,1,0.5)) # https://docs.makie.org/stable/explanations/inspector/index.html
    display(GLMakie.Screen(), fig) # note the window only lasts as long as the julia session
  end
  # dev1129()

  function devInvoluteXYs()
    # I need to thoroughly and graphically test functions relating involute angle to xy position to determine what broke when inserting UnitTypes
    # First, calcThetaHeight is not finding the correct angle, possibly due to a confusion on whether the radial height is absolute or relative to rb.

    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    # (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)
    un = Inch

    fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
    axs = Axis(fig[1,1], xlabel="X[$(un(1).unit)]", ylabel="Y[$(un(1).unit)]", aspect=DataAspect())

    alBase = Radian(1)
    gm = Radian(1.2)
    Gears.InvoluteTooth.plotArc(axs=axs, r=g.base.measure/2, a0=Radian(0), a1=Radian(2), label="base", linecolor=:black)
    Gears.InvoluteTooth.plotArc(axs=axs, r=g.root.measure/2, a0=Radian(0), a1=Radian(2), label="root", linecolor=:blue)
    Gears.InvoluteTooth.plotArc(axs=axs, r=g.pitch.measure/2, a0=Radian(0), a1=Radian(2), label="pitch", linecolor=:green)
    Gears.InvoluteTooth.plotArc(axs=axs, r=g.outside.measure/2, a0=Radian(0), a1=Radian(2), label="outside", linecolor=:red)

    Gears.InvoluteTooth.plotInvolute(axs=axs, bd=g.base, al=alBase, thMax=Radian(2))

    Gears.InvoluteTooth.plotInvoluteConstruction(axs=axs, bd=g.base, al=alBase, th= Gears.InvoluteTooth.findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=(g.pitch/2).measure ))
    Gears.InvoluteTooth.plotInvoluteConstruction(axs=axs, bd=g.base, al=alBase, th= Gears.InvoluteTooth.findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=(g.outside/2).measure ))


    DataInspector(fig, transparency=true, backgroundcolor=RGBAf(1,1,1,0.5)) # https://docs.makie.org/stable/explanations/inspector/index.html
    display(GLMakie.Screen(), fig) # note the window only lasts as long as the julia session

  end
  devInvoluteXYs()

end
;