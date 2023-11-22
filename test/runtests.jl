using TestItemRunner
@run_package_tests 


# using Test
# using Plots
# gr()
# # plotly()
# # pythonplot()
# using Gears
# using UnitTypes



# @testset "functions" begin

#   @test Gears.pitchDiameter(30, Gears.Pitch(24)) ≈ Diameter(Inch(1.250))
#   @test Gears.pitchDiameter(Gears.Pitch(24), 30) ≈ Diameter(Inch(1.250))

#   n = 30
#   pa = Degree(20)
#   dp = Gears.Pitch(24)
#   pd = Gears.pitchDiameter(n, dp)



#   bd = Gears.baseDiameter(pd, pa)
#   @test isapprox(bd, Inch(1.1746), atol=1e-3)

#   dp = Gears.diametralPitch(n, pd)
#   @test isapprox(dp, Gears.Pitch(24), atol=1e-3)

#   ad = Gears.addendum(dp)
#   @test isapprox(ad, Inch(1/24), atol=1e-3)

#   dd = Gears.dedendum(dp)
#   @test isapprox(dd, Inch(1.25/24), atol=1e-3)

#   od = Gears.outsideDiameter(pd, ad) 
#   @test isapprox(od, Diameter(Inch(1.3333)), atol=1e-3)

#   rd = Gears.rootDiameter(pd, dd)
#   @test isapprox(rd, pd-Diameter(2*dd), atol=1e-3)

#   bi = Gears.baseInterval(n, bd)
#   @test typeof(bi) <: UnitTypes.AbstractExtent
#   @test isapprox(bi, π*bd.value/n)

# end

# # include("test_GearANSI.jl")

# # include("test_InvoluteTooth.jl")
