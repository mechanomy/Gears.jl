module Dev
  using UnitTypes
  using Gears
  using GLMakie

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
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    # (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)

    


  end
  dev1129()

end
;