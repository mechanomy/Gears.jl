using Test
using Plots
using Gears
using UnitTypes


@testset "functions" begin
  n = 30
  pa = Degree(20)
  pd = Diameter(Inch(1.250))

  bd = Gears.baseDiameter(pd, pa)
  @test isapprox(bd, Inch(1.1746), atol=1e-3)

  dp = Gears.diametralPitch(n, pd)
  @test isapprox(dp, Gears.Pitch(24), atol=1e-3)

  ad = Gears.addendum(dp)
  @test isapprox(ad, Inch(1/24), atol=1e-3)

  dd = Gears.dedendum(dp)
  @test isapprox(dd, Inch(1.25/24), atol=1e-3)

  od = Gears.outsideDiameter(pd, ad) 
  @test isapprox(od, Diameter(Inch(1.3333)), atol=1e-3)

  rd = Gears.rootDiameter(pd, dd)
  @test isapprox(rd, pd-Diameter(2*dd), atol=1e-3)

  bi = Gears.baseInterval(n, bd)
  @test typeof(bi) <: UnitTypes.AbstractExtent
  @test isapprox(bi, π*bd.value/n)

end

@testset "GearANSI" begin

  # https://shop.sdp-si.com/catalog/product/?id=S10C9Z-024H030
  nTeeth = 30
  pressureAngleDegree = Degree(20)
  # diametralPitch = 24
  pitchDiameterInch = Diameter(Inch(1.250))
  boreInch = Inch(0.4998)
  # hubDiameterInch = 0.75
  # faceWidthInch = 0.3750
  # agmaQualityClass = 10

  # @show convert(Float64, boreInch)
  # @show convert(Float64, pressureAngleDegree)
  # @show convert(Float64, pitchDiameterInch)

  @show s = Gears.GearANSI( n=nTeeth, pitch=pitchDiameterInch, pressure=pressureAngleDegree )

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
  @test true

  # make an example for plotting and math
  s = Gears.GearANSI( n=30, pitch=Diameter(Inch(1.250)), pressure=Degree(20) )
  rod = convert(Float64,s.outside.value/2) # make a strip()?
  rpd = convert(Float64,s.pitch.value/2)
  rbd = convert(Float64,s.base.value/2)
  rrd = convert(Float64,s.root.value/2)

  p = plot()
  p = Gears.InvoluteTooth.plotBaseSection( r=rod, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)
  p = Gears.InvoluteTooth.plotBaseSection( r=rpd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:black, linestyle=:dash)
  p = Gears.InvoluteTooth.plotBaseSection( r=rbd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)
  p = Gears.InvoluteTooth.plotBaseSection( r=rrd, gm=0, al=0, bt=π/4, thExtra=0.5, n=100, p=p, linecolor=:gray, linestyle=:dash)

  i = 3
  gm=2*π/s.nTeeth*(i-1)
  # al = gm - (π/2/s.nTeeth + tan(s.pressure)-s.pressure)  ## need iterate
  al = gm - (π/2/s.nTeeth + tan(convert(Radian,s.pressure).value)-convert(Radian,s.pressure).value )
  ht = rod-rrd
  althPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=al, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=al, thMax=althPdOd, p=p, linecolor=:blue, linestyle=:solid)
  p = Gears.InvoluteTooth.plotInvoluteConstruction( r=rrd, al=al, th=althPdOd, p=p, linecolor=:gray, linestyle=:dash)
  # p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.6, aMax=, aMin=al, label="γ", linecolor=:gray, linestyle=:solid, fontsize=10)
  p = Gears.InvoluteTooth.plotGamma( p=p, gm=gm, r=rrd, ht=ht, linecolor=:green, linestyle=:solid )
  p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.4, aMax=gm, aMin=0, label="γ", linecolor=:gray, linestyle=:solid, fontsize=10)
  p = Gears.InvoluteTooth.plotAlpha( p=p, r=rrd, al=al, linecolor=:gray, linestyle=:solid)
  p = Gears.InvoluteTooth.drawAngleArc( p=p, r=rrd*.5, aMax=al, aMin=0, label="α", linecolor=:gray, linestyle=:solid, fontsize=10)

  # bt = gm + (π/2/s.nTeeth + tan(s.pressure)-s.pressure)
  bt = gm + (π/2/s.nTeeth + tan(convert(Radian,s.pressure).value)-convert(Radian,s.pressure).value )
  btthPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=bt, ht=ht )
  p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=bt, thMax=btthPdOd, p=p, linecolor=:red, linestyle=:solid)

  # for i in 1:s.nTeeth
  #   gm=2*π/s.nTeeth*(i-1)
  #   al = gm - (π/2/s.nTeeth + tan(s.pressure)-s.pressure) 
  #   althPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=al, ht=ht )
  #   p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=al, thMax=althPdOd, p=p, linecolor=:blue, linestyle=:solid)

  #   bt = gm + (π/2/s.nTeeth + tan(s.pressure)-s.pressure) 
  #   btthPdOd = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=bt, ht=ht )
  #   p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=bt, thMax=btthPdOd, p=p, linecolor=:red, linestyle=:solid)
  # end
  display(p)
  # # savefig("julia/Gears/on230630.svg")

end;
