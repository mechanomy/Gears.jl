using Test
using UnitTypes


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
  # p = plot()
  # p = plot(p,xs, ys, linecolor=:blue, linestyle=:dash, linesize=1, markersize=3, markercolor=:red, aspect_ratio=:equal, linealpha=0.5)
  # display(p)
end



# @testset "InvoluteTooth" begin
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=0.3, ht=20 ), 1.238, atol=1e-3)
#   @test isapprox( Gears.InvoluteTooth.calcThetaHeight(r=96, gm=deg2rad(56), al=1.27, ht=20 ), 0.550, atol=1e-3)
# end


# 1989 hindhede lists properties of involutes that might be tested

