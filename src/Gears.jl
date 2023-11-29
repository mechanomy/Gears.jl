module Gears
  # References:
  # 2004_oberg ~p1599
  # 1989_Hindhede_SpurGearFundamentals - has good diagrams
  # 2012_Dooner_KinematicGeometryOfGearing ch2 - relates tooth width to involute start angle

  using UnitTypes
  using TestItems
  using DocStringExtensions

  #inverse length not part of UnitTypes yet
  @makeBaseMeasure InverseLength PerMeter "m^-1"
  @makeDerivedMeasure PerInch "in^-1" 1/0.0254 PerMeter # 39.3700 in/meter
  Base.:/(x::T, y::U) where {T<:Number, U<:AbstractLength} = PerInch(x / convert(Inch, y).value )
  Base.:/(x::T, y::U) where {T<:Number, U<:AbstractInverseLength} = Inch(x / convert(PerInch, y).value )
  # Base.:*(x::T, y::U) where {T<:AbstractLength, U<:AbstractInverseLength} = x.value*x.toBase * y.value*y.toBase
  @testitem "PerInch" begin
    using UnitTypes
    @test PerMeter(2).value ≈ 2
    @test PerInch(1) ≈ PerMeter(1/.0254) # 1in/0.0254m
    @test 1/Inch(2) ≈ PerInch(0.5)
    @test isapprox( 1/Inch(3) , PerInch(1/3), atol=1e-3)
    @test 1/PerInch(2) ≈ Inch(0.5)

    # @test Inch(3) * PerInch(3) ≈ 1 # 3/inch*3inch 
  end

  # Gear-specific types:
  """
  the diameter of the pitch circle
  """
  #not documentable
  @makeDimension PitchDiameter Inch # diameter of the pitch circle
  Base.isapprox(x::T, y::U; atol::Real=0, rtol::Real=atol) where {T<:PitchDiameter, U<:PitchDiameter} = isapprox(x.measure, convert(typeof(x.measure), y.measure), atol=atol, rtol=rtol) # convert both into the same unit
  @testitem "PitchDiameter" begin
    using UnitTypes
    @test isapprox(PitchDiameter(Millimeter(86.36)), PitchDiameter(Inch(3.4)), atol=1e-3 )

    nt = 10
    pd = PitchDiameter(Inch(3.4))
    dp = DiametralPitch(nt,pd)
    pd2 = PitchDiameter(nt, dp)
    @test isapprox(pd, pd2, atol=1e-3)
  end

  @makeDimension DiametralPitch PerInch # = Nteeth/PitchDiameter, the number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch
  DiametralPitch(nTeeth::Int, pd::PitchDiameter) = DiametralPitch( nTeeth / convert(Inch, pd.measure)) # macro defines DiametralPitch{PerInch}(num), but this is the actual definition
  # instead of using the macro, define by hand so that I can have a custom constructor
  # abstract type AbstractDiametralPitch <: AbstractDimension end
  """
  the number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch
  """
  # struct DiametralPitch{T<:AbstractInverseLength} <: AbstractDiametralPitch
  #   measure::T
  #   DiametralPitch(nTeeth::T, pd::U) where {T<:Int, U<:PitchDiameter} = new{PerInch}( nTeeth / convert(Inch,pd.measure) )
  # end
  # export DiametralPitch
  # Base.:/(x::T, y::U) where {T<:Number, U<:AbstractDiametralPitch} = PitchDiameter(x / y.measure) # not true as it does not know nTeeth
  @testitem "DiametralPitch" begin
    using UnitTypes
    dp = DiametralPitch(30, PitchDiameter(Inch(1.250)) )
    @test dp.measure ≈ PerInch(24)

    dp = DiametralPitch(10, PitchDiameter(Millimeter(50.8)) )
    @test dp.measure ≈ PerInch(5)
  end

  #needs to follow DP definition
  PitchDiameter(dp::DiametralPitch, nTeeth::Int) = PitchDiameter( nTeeth / dp.measure )
  PitchDiameter(nTeeth::Int, dp::DiametralPitch) = PitchDiameter( nTeeth / dp.measure ) # positions shouldn't matter!

  @makeDimension Addendum Inch
  """
  height of tooth above pitch circle
  """
  Addendum(dp::DiametralPitch) = Addendum( 1/dp.measure )

  @makeDimension Dedendum Inch
  """
  depth of tooth below pitch circle
  """
  Dedendum(dp::DiametralPitch) = Dedendum( 1.250/dp.measure )

  @testitem "Addendum&Dedendum" begin
    using UnitTypes
    dp = DiametralPitch(10, PitchDiameter(Inch(2)))
    @test Addendum(dp) ≈ Addendum(Inch(0.2))
    @test Dedendum(dp) ≈ Dedendum(Inch(0.25))
  end

  # outsideDiameter(pd::PitchDiameter, addendum) = pd + 2*addendum
  # vs making a type and a custom constructor?
  @makeDimension OutsideDiameter Inch
  """
  Outside diameter is the maximum diameter of the gear teeth
  """
  OutsideDiameter(pd::PitchDiameter, ad::Addendum) = OutsideDiameter( (pd + 2*ad).measure ) # () is a PitchDiameter, so we need to get the Inch measure

  @testitem "OutsideDiameter" begin
    using UnitTypes
    nt = 30
    pd = PitchDiameter(Inch(1.2500))
    od = OutsideDiameter(Inch(1.333))
    dp = DiametralPitch(nt, pd)
    ad = Addendum(dp)
    @test isapprox(od.measure, OutsideDiameter(pd, ad).measure, atol=1e-3)
    @test isapprox(od - pd, Inch(1/3-1/4), atol=1e-3)
    @test isapprox(OutsideDiameter(Inch(1))-BaseDiameter(Millimeter(25.4)), Meter(0), atol=1e-3 ) # just check that we still handle unit conversions correctly
  end


  @makeDimension BaseDiameter Inch
  """
  tooth involutes are perpendicular to this circle
  in meshing gears, connecting two base circles with a tangent line will cross the involute teeth at the point of contact
  """
  BaseDiameter(pd::PitchDiameter, pressureAngle::U where U<: AbstractAngle) = BaseDiameter( (pd * cos(pressureAngle)).measure)
  @testitem "BaseDiameter" begin
    using UnitTypes
    pd = PitchDiameter(Inch(1.2500))
    pa = Degree(20)
    bd = BaseDiameter(pd, pa)
    @test isapprox(bd.measure, Inch(1.1746), atol=1e-3)
  end

  @makeDimension RootDiameter Inch
  """
  Diameter of the circle touching the bottom lands between teeth
  """
  RootDiameter(pd::PitchDiameter, dd::Dedendum) = RootDiameter((pd - 2*dd).measure)
  @testitem "RootDiameter" begin
    using UnitTypes
    nt = 30
    pd = PitchDiameter(Inch(1.2500))
    dp = DiametralPitch(nt, pd)
    dd = Dedendum(dp)
    rd = RootDiameter(pd, dd)
    @test isapprox(rd.measure, Inch(1.14583), atol=1e-3)
  end



  """
  The number of teeth on the gear described by `dp` and `pd`; note that this does not ensure that `dp` and `pd` are consistent, rounding to the nearest tooth!
  """
  nTeeth(dp::DiametralPitch, pd::PitchDiameter) = Int64(round( convert(PerInch,dp.measure).value * convert(Inch, pd.measure).value ))
  # nTeeth(pitchD::T, circularPitch) = π*pitchDiameter/circularPitch
  export nTeeth
  @testitem "nTeeth" begin
    using UnitTypes
    nt = 30
    pd = PitchDiameter(Inch(1.2500))
    dp = DiametralPitch(nt, pd)
    @test nTeeth(dp, pd) ≈ nt
  end


  # """
  # distance between successive involutes measured along the segment of the base circle, mating teeth must have the same base pitch length
  # """
  # # basePitch(nTeeth::Int, baseDiameter::Diameter) = π*baseDiameter/nTeeth # length/tooth
  # baseInterval(nTeeth::Int, baseDiameter::Diameter) = π*convert(typeof(baseDiameter).parameters[1],baseDiameter)/nTeeth # length/tooth
  

  @makeDimension CircularPitch Inch
  """
    the arc distance between similar points on the adjacent teeth in linear units along the arc segment
  """
  CircularPitch(pd::PitchDiameter, nTeeth) = π*pd / nTeeth
  # circularPitch(diametralPitch) = π/diametralPitch
  @testitem "CircularPitch" begin
    using UnitTypes
    pd = PitchDiameter(Inch(1.250))
    nt = 30
    cp = CircularPitch(pd, nt)
    @test isapprox( cp.measure, Inch(0.130899), atol=1e-3 )
  end

  @makeDimension GearModule Millimeter
  """
  module is the amount of pitch diameter per tooth, higher module larger tooth
  """
  GearModule(pd::PitchDiameter, nTeeth) = GearModule( convert(Millimeter, pd.measure) / nTeeth )
  @testitem "GearModule" begin
    using UnitTypes
    # S10T05M150S0308
    mo = GearModule(Millimeter(0.5))
    nt = 150
    pa = Degree(20)
    pd = PitchDiameter(Millimeter(75.00))
    @test GearModule(pd, nt) ≈ mo
  end


  abstract type AbstractGear end

  include("GearANSI.jl")
  include("InvoluteTooth.jl")

end #Gears

