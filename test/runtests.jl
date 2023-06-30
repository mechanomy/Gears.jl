using Pkg
Pkg.activate( normpath(joinpath(@__DIR__, "..")) ) #activate this package

using Test
using Plots
using Gears
using UnitTypes


@testset "SpurANSI" begin

  # https://shop.sdp-si.com/catalog/product/?id=S10C9Z-024H030
  nTeeth = 30
  pressureAngleDegree = 20
  diametralPitch = 24
  # pitchDiameterInch = 1.250
  # boreInch = 0.4998
  # hubDiameterInch = 0.75
  # faceWidthInch = 0.3750
  # agmaQualityClass = 10

  @show s = Gears.SpurANSI( nTeeth, Degree(20), Gears.Pitch(24))

end


# @testset "InvoluteTooth" begin
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=0.3, ht=20 ), 1.238, atol=1e-3)
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=1.27, ht=20 ), 0.550, atol=1e-3)
# end


# @testset "InvoluteTooth calc" begin
#   p = Gears.InvoluteTooth.plotTooth( r=96, gm=deg2rad(45), dl=deg2rad(56), ht=20 )
#   # savefig("on220826_involutes_toothHeight.svg")
#   display(p)
#   @test true
# end

#1989 hindhede lists properties of involutes that might be tested

@testset "gear info to involute shape" begin

  # make an example for plotting and math
  nTeeth = 30
  pressureAngleDegree = 20
  diametralPitch = 24
  pitchDiameter = nTeeth/diametralPitch
  rpd = pitchDiameter/2

  addendum = 1/diametralPitch
  dedendum = 1.250/diametralPitch

  rootDiameter = pitchDiameter - 2*dedendum
  rrd = rootDiameter/2

  baseDiameter = pitchDiameter*cos(deg2rad(pressureAngleDegree))
  rbd = baseDiameter/2

  outsideDiameter = pitchDiameter + 2*addendum
  rod = outsideDiameter/2

  # gm = deg2rad(45) # tooth symmetry angle
  basePitch = π*baseDiameter/nTeeth # = distance between successive involutes
  basePitchAngle = 2*π/nTeeth # angle of one tooth+groove

  p = plot()
  p = Gears.InvoluteTooth.plotBaseSection( r=rod, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)
  p = Gears.InvoluteTooth.plotBaseSection( r=rpd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:black, linestyle=:dash)
  # p = Gears.InvoluteTooth.plotBaseSection( r=rbd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)
  p = Gears.InvoluteTooth.plotBaseSection( r=rrd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)

  i = 3
  @show gm=2*π/nTeeth*(i-1)
  # dl = thicknessAtPitchDiameter/(pitchDiameter/2)
  # al = gm - dl/2
  @show al = gm - (π/2/nTeeth + tan(deg2rad(pressureAngleDegree))-deg2rad(pressureAngleDegree) )
  ht = rod-rrd
  althPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=al, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=al, thMax=althPdOd, p=p, linecolor=:blue, linestyle=:solid)
  p = Gears.InvoluteTooth.plotInvoluteConstruction( r=rrd, al=al, th=althPdOd, p=p, linecolor=:gray, linestyle=:dash)
  # p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.6, aMax=, aMin=al, label="γ", linecolor=:gray, linestyle=:solid, fontsize=10)
  p = Gears.InvoluteTooth.plotGamma( p=p, gm=gm, r=rrd, ht=ht, linecolor=:green, linestyle=:solid )
  p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.4, aMax=gm, aMin=0, label="γ", linecolor=:gray, linestyle=:solid, fontsize=10)
  p = Gears.InvoluteTooth.plotAlpha( p=p, r=rrd, al=al, linecolor=:gray, linestyle=:solid)
  p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.5, aMax=al, aMin=0, label="α", linecolor=:gray, linestyle=:solid, fontsize=10)
  

  @show bt = gm + (π/2/nTeeth + tan(deg2rad(pressureAngleDegree))-deg2rad(pressureAngleDegree) )
  btthPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=bt, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=bt, thMax=btthPdOd, p=p, linecolor=:red, linestyle=:solid)

  # i = 4
  for i in 1:nTeeth
  @show gm=2*π/nTeeth*(i-1)
  @show al = gm - (π/2/nTeeth + tan(deg2rad(pressureAngleDegree))-deg2rad(pressureAngleDegree) )
  althPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=al, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=al, thMax=althPdOd, p=p, linecolor=:blue, linestyle=:solid)

  @show bt = gm + (π/2/nTeeth + tan(deg2rad(pressureAngleDegree))-deg2rad(pressureAngleDegree) )
  btthPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=bt, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=bt, thMax=btthPdOd, p=p, linecolor=:red, linestyle=:solid)

  # p = plot()
  # for i in 1:nTeeth
  #   p = Gears.InvoluteTooth.plotTooth(r=rootDiameter/2, gm=2*π/nTeeth*(i-1), dl=dl, ht=(outsideDiameter-rootDiameter)/2, p=p)
  end
  display(p)
  # savefig("julia/Gears/on230630.svg")

end
