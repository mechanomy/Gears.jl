using Test
using Plots
gr()
# plotly()
# pythonplot()
using Gears
using UnitTypes



@testset "functions" begin
  n = 30
  pa = Degree(20)
  dp = Gears.Pitch(24)

  pd = Gears.pitchDiameter(n, dp) #Diameter(Inch(1.250))
  pd = Gears.pitchDiameter(dp, n) #Diameter(Inch(1.250))

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
  diametralPitch = Gears.Pitch(24)
  pitchDiameterInch = Diameter(Inch(1.250))
  addendum = Inch(1/24)
  dedendum = Inch(1.25/24)
  outside = Diameter(Inch(1.333))
  base = Diameter(Inch(1.1746))
  root = Diameter(Inch(1.1458))
  baseInterval = Inch(0.12300)
  # boreInch = Inch(0.4998)
  # hubDiameterInch = 0.75
  # faceWidthInch = 0.3750
  # agmaQualityClass = 10


  #constructor permutations
  gref = Gears.GearANSI( nTeeth, pitchDiameterInch, pressureAngleDegree, diametralPitch, addendum, dedendum, outside, base, root, baseInterval )
  @test typeof(gref) <: Gears.GearANSI

  @test typeof(Gears.GearANSI(nTeeth, diametralPitch)) <: Gears.GearANSI
  @test typeof(Gears.GearANSI(nTeeth, diametralPitch, pressureAngleDegree)) <: Gears.GearANSI
  @test typeof(Gears.GearANSI(diametralPitch, nTeeth, pressureAngleDegree)) <: Gears.GearANSI
  @test typeof(Gears.GearANSI( 30, Diameter(Inch(1.250)), Degree(20) )) <: Gears.GearANSI
  # @test typeof(Gears.GearANSI( 30, Diameter(MilliMeter(31.75)), Degree(20) )) <: Gears.GearANSI
end

@testset "calcThetaHeight" begin
  g = Gears.GearANSI( 30, Diameter(Inch(1.250)), Degree(20) )

  # @show thetaH = Gears.InvoluteTooth.calcThetaHeight(r=convert(Radius,g.base), gm=Degree(45), al=Degree(40), toothHeight=convert(Meter,(g.outside-g.base)/2) )
  # @test isapprox(thetaH, 7.9524, atol=1e-3)
  # @test typeof(thetaH) <: AbstractExtent
end

@testset "gear to plane points" begin
  # g = Gears.GearANSI( n=30, pitch=Diameter(Inch(1.250)), pressure=Degree(20) ) #sdpsi_S10C9Z-024H030
  g = Gears.GearANSI( 70, Diameter(Inch(2.9167)), Degree(20) ) #sdpsi_S1268Z-024A070
  # @show g = Gears.GearANSI( Gears.Pitch(24), Diameter(Inch(100*2/25.4)), Degree(20) ) not working, need to do coercion
  # @show g = Gears.GearANSI(Gears.Pitch(24), 130, Degree(20) )
  # @show g = Gears.GearANSI(Gears.Pitch(24), 3000, Degree(20) )
  # g = Gears.GearANSI(Gears.Pitch(24), 189, Degree(20) ) # 100mm radius = 7.87in diameter = 188.8 teeth=> 189

  nPerTooth = 50
  nPerInvolute = convert(Int64,round(nPerTooth/2)) #this may be better defined by physical spacing...
  Gears.InvoluteTooth.writeToothProfilePoints(g, nPerTooth=20)
  (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)
  p = plot()
  p = plot(p,xs, ys, linecolor=:blue, linestyle=:dash, linesize=1, markersize=3, markercolor=:red, aspect_ratio=:equal, linealpha=0.5)
  display(p)
end



# @testset "InvoluteTooth" begin
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=0.3, ht=20 ), 1.238, atol=1e-3)
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=1.27, ht=20 ), 0.550, atol=1e-3)
# end


# 1989 hindhede lists properties of involutes that might be tested

