# struct and functions for working with ANSI-spec sprug gears.
# As the full specification is not available, included here are the most basic and common parameters and calculations.

"""
  Struct to represent ANSI-specification spur gears.
"""
struct GearANSI <: AbstractGear
  nTeeth::Int
  pitch::PitchDiameter #again, I don't care what the unit of measure is just as long as it's convertible with a length. I am a little surprised that specification is not required...
  pressure::U where U<:UnitTypes.AbstractAngle
  diametral::DiametralPitch
  addendum::Addendum
  dedendum::Dedendum
  outside::OutsideDiameter
  base::BaseDiameter
  root::RootDiameter

  # baseInterval::A where A<:UnitTypes.AbstractLength
  # bore::UnitTypes.AbstractDiameter
  # circularPitch::Pitch
  # wholeDepth::UnitTypes.AbstractLength
  # workingDepth::UnitTypes.AbstractLength
  # clearance::UnitTypes.AbstractLength
end

function GearANSI(dp::DiametralPitch, nTeeth::Int, pa::UnitTypes.Degree=UnitTypes.Degree(20)) 
  pd = PitchDiameter(dp, nTeeth)
  bd = BaseDiameter(pd,pa)
  ad = Addendum(dp)
  dd = Dedendum(dp)
  od = OutsideDiameter(pd, ad)
  rd = RootDiameter(pd, dd)
  return GearANSI( nTeeth, pd, pa, dp, ad, dd, od, bd, rd)
end

function GearANSI(pd::PitchDiameter, nTeeth::Int, pa::UnitTypes.Degree=UnitTypes.Degree(20)) 
  dp = DiametralPitch(nTeeth, pd)
  bd = BaseDiameter(pd,pa)
  ad = Addendum(dp)
  dd = Dedendum(dp)
  od = OutsideDiameter(pd, ad)
  rd = RootDiameter(pd, dd)
  return GearANSI( nTeeth, pd, pa, dp, ad, dd, od, bd, rd)
end
export GearANSI

@testitem "GearANSI" begin
  using UnitTypes
  # https://shop.sdp-si.com/catalog/product/?id=S10C9Z-024H030
  nTeeth = 30
  pa = Degree(20)
  dp = DiametralPitch(PerInch(24))
  pd = PitchDiameter(Inch(1.250))
  ad = Addendum(Inch(1/24))
  dd = Dedendum(Inch(1.25/24))
  od = OutsideDiameter(Inch(1.333))
  bd = BaseDiameter(Inch(1.1746))
  rd = RootDiameter(Inch(1.1458))
  # baseInterval = Inch(0.12300)
  # boreInch = Inch(0.4998)
  # hubDiameterInch = 0.75
  # faceWidthInch = 0.3750
  # agmaQualityClass = 10

  #constructor permutations
  @test typeof(GearANSI( nTeeth, pd, pa, dp, ad, dd, od, bd, rd)) <: GearANSI
  @test typeof(GearANSI(dp, nTeeth, pa)) <: GearANSI
  @test typeof(GearANSI(dp, nTeeth)) <: GearANSI
  @test typeof(GearANSI(pd, nTeeth)) <: GearANSI
end

"""
    pulley2String(p::PlainPulley) :: String
  Returns a descriptive string of the given PlainPulley `p` of the form:
    PlainPulley[struct] @ [1.000mm,2.000mm] r[3.000mm] arrive[57.296°] depart[114.592°] aWrap[57.296°] lWrap[3.000mm]"
"""
function gear2String(g::GearANSI)::String 
  # return @sprintf("GearANSI:\n nTeeth: %d \n pressureAngleDegree: %3.3f[°] \n pitchDiameterInch: %3.3f[in], diametralPitch: %3.3f[1/in]",
  #   g.nTeeth,
  #   convert(Float64, g.pressureAngle.value),
  #   convert(Float64, g.pitchDiameter.value),
  #   convert(Float64, g.diametralPitch.value)
  # )
  return "GearANSI: $(g.nTeeth)tooth $(g.pitch)"
end
@testitem "gear2String" begin
  using UnitTypes
  g = GearANSI(PitchDiameter(Inch(1.2500)), 30, Degree(20))
  # @test Gears.gear2String(g) === "GearANSI: 30tooth Gears.PitchDiameter(1.25in)"
  @test Gears.gear2String(g) === "GearANSI: 30tooth Gears.GearDimensions.PitchDiameter(1.25in)"
end

"""
  Function to `show()` a GearANSI via [`gear2String()`](@ref)
"""
function Base.show(io::IO, g::GearANSI)
  print(io, gear2String(g))
end


