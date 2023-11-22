module Dev
  using UnitTypes
  using Gears
  using GLMakie

  GLMakie.closeall()
  GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config


  # g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
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

  # Gears.InvoluteTooth.plotBaseSection(axs=axs, r=r, al=al, bt=bt)
  # Gears.InvoluteTooth.plotGamma(axs=axs, gm=gm, r=r, ht=0.5 )
  # Gears.InvoluteTooth.plotAlpha(axs=axs, r=r, al=al)
  # Gears.InvoluteTooth.plotInvolute(axs=axs, r=r, al=al, thMax=0.2)
  # Gears.InvoluteTooth.plotInvoluteConstruction(axs=axs, r=r, al=al, th=1.2)
  # Gears.InvoluteTooth.drawAngleArc(axs=axs, r=r, aMax=1.5, aMin=al)
  # Gears.InvoluteTooth.drawToothTop(axs=axs, r=r, al=al, thal=0.2, bt=bt, thbt=0.2)

  Gears.InvoluteTooth.plotTooth(axs=axs, r=r, gm=gm, dl=bt-al, ht=ht)

  display(GLMakie.Screen(), fig)
  # display(fig)

end