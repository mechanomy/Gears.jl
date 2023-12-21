module GearDimensions

  #references:
  # https://en.wikipedia.org/wiki/List_of_gear_nomenclature
  # https://khkgears.net/new/gear_knowledge/gear_technical_reference/calculation_gear_dimensions.html


  using TestItems
  using UnitTypes

  @makeDimension PitchDiameter Inch # diameter of the pitch circle

  @testitem "PitchDiameter" begin
    using UnitTypes
    @test isapprox(PitchDiameter(MilliMeter(86.36)), PitchDiameter(Inch(3.4)), atol=1e-3 )
  end

  @makeBaseMeasure InverseLength PerMeter "m^-1"  #inverse length not part of UnitTypes yet
  # @makeDerivedMeasure PerInch "in^-1" 1/0.0254 PerMeter # 39.3700 in/meter
  @deriveMeasure PerMeter(1/0.0254) = PerInch(1) "in^-1"
  Base.:/(x::T, y::U) where {T<:Real, U<:AbstractLength} = PerInch(x / convert(Inch, y).value )
  Base.:/(x::T, y::U) where {T<:Real, U<:AbstractInverseLength} = Inch(x / convert(PerInch, y).value )

  Base.round(x::T, r::RoundingMode=RoundNearest; digits::Integer=0, base=10) where T<:AbstractMeasure = T( round(x.value, r, digits=digits, base=base) )
  @testitem "PerInch" begin
    using UnitTypes
    @test PerMeter(2).value ≈ 2
    @test PerInch(1) ≈ PerMeter(1/.0254) # 1in/0.0254m, so 10qty/inch = 393.7qty/meter
    @test 1/Inch(2) ≈ PerInch(0.5)
    @test isapprox( 1/Inch(3) , PerInch(1/3), atol=1e-3)
    @test 1/PerInch(2) ≈ Inch(0.5)
    # @test Inch(3) * PerInch(3) ≈ 1 # 3/inch*3inch 
  end

  # Since DiametralPitch is the number of teeth in an inch, it is defined to be an integer
  @makeDimension DiametralPitch PerInch # = Nteeth/PitchDiameter, the number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch
  """
    $TYPEDSIGNATURES

    The number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch.
  """
  DiametralPitch(nTeeth::Int, pd::PitchDiameter) = DiametralPitch( round(nTeeth / convert(Inch, pd.measure))) # macro defines DiametralPitch{PerInch}(num), but this is the actual definition
  @testitem "DiametralPitch" begin
    using UnitTypes
    dp = DiametralPitch(30, PitchDiameter(Inch(1.250)) )
    @test dp.measure ≈ PerInch(24)

    @show dp = DiametralPitch(30, PitchDiameter(Inch(1.23)) )
    @test dp.measure ≈ PerInch(24)

    dp = DiametralPitch(10, PitchDiameter(MilliMeter(50.8)) )
    @test dp.measure ≈ PerInch(5)
  end

  #needs to follow DP definition
  """
    $TYPEDSIGNATURES

    The diameter of the pitch circle which, when engaged with another gear has the same angular velocity as that other gear's pitch circle
  """
  PitchDiameter(dp::DiametralPitch, nTeeth::Int) = PitchDiameter( nTeeth / dp.measure )
  PitchDiameter(nTeeth::Int, dp::DiametralPitch) = PitchDiameter( nTeeth / dp.measure ) # positions shouldn't matter!

  @makeDimension Addendum Inch
  """
    $TYPEDSIGNATURES
    Height of tooth above pitch circle.
  """
  Addendum(dp::DiametralPitch) = Addendum( 1/dp.measure )

  @makeDimension Dedendum Inch
  """
    $TYPEDSIGNATURES

    Depth of tooth below pitch circle.
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
    $TYPEDSIGNATURES

    Outside diameter is the maximum diameter of the gear teeth.
  """
  OutsideDiameter(pd::PitchDiameter, ad::Addendum) = OutsideDiameter(pd + 2*ad) # () is a PitchDiameter, so we need to get the Inch measure
  Base.:-(x::OutsideDiameter, y::PitchDiameter) = x.measure-y.measure
  PitchDiameter(od::OutsideDiameter, nTeeth::Int) = od / (1 + 2/nTeeth) # coerce to int dp?

  @testitem "OutsideDiameter" begin
    using UnitTypes
    nt = 30
    pd = PitchDiameter(Inch(1.2500))
    od = OutsideDiameter(Inch(1.333))
    dp = DiametralPitch(nt, pd)
    ad = Addendum(dp)
    @test isapprox(od.measure, OutsideDiameter(pd, ad).measure, atol=1e-3)
    @test isapprox(od - pd, Inch(1/3-1/4), atol=1e-3)
    @test isapprox(OutsideDiameter(Inch(1))-BaseDiameter(MilliMeter(25.4)), Meter(0), atol=1e-3 ) # just check that we still handle unit conversions correctly
  end

  @makeDimension BaseDiameter Inch
  """
    $TYPEDSIGNATURES

    Tooth involutes are perpendicular to this circle.
    In meshing gears, connecting two base circles with a tangent line will cross the involute teeth at the point of contact.
  """
  BaseDiameter(pd::PitchDiameter, pressureAngle::U where U<: AbstractAngle) = BaseDiameter( pd.measure * cos(pressureAngle))
  Base.:-(x::OutsideDiameter, y::BaseDiameter) = x.measure-y.measure

  @testitem "BaseDiameter" begin
    using UnitTypes
    pd = PitchDiameter(Inch(1.2500))
    pa = Degree(20)
    bd = BaseDiameter(pd, pa)
    @test isapprox(bd.measure, Inch(1.1746), atol=1e-3)
  end

  @makeDimension RootDiameter Inch
  """
    $TYPEDSIGNATURES

    Diameter of the circle touching the bottom lands between teeth.
  """
  RootDiameter(pd::PitchDiameter, dd::Dedendum) = RootDiameter(pd.measure - 2*dd)
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
    $TYPEDSIGNATURES

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


  @makeDimension CircularPitch Inch
  """
    $TYPEDSIGNATURES
    The arc distance between similar points on the adjacent teeth in linear units along the arc segment.
  """
  CircularPitch(pd::PitchDiameter, nTeeth) = CircularPitch(Inch(π*pd.measure.value / nTeeth))
  # circularPitch(diametralPitch) = π/diametralPitch
  @testitem "CircularPitch" begin
    using UnitTypes
    pd = PitchDiameter(Inch(1.250))
    nt = 30
    cp = CircularPitch(pd, nt)
    @test isapprox( cp.measure, Inch(0.130899), atol=1e-3 )
  end

  @makeDimension GearModule MilliMeter
  """
    $TYPEDSIGNATURES

    Module is the amount of pitch diameter per tooth, higher module larger tooth.
  """
  GearModule(pd::PitchDiameter, nTeeth) = GearModule( convert(MilliMeter, pd.measure) / nTeeth )
  @testitem "GearModule" begin
    using UnitTypes
    # S10T05M150S0308
    mo = GearModule(MilliMeter(0.5))
    nt = 150
    pa = Degree(20)
    pd = PitchDiameter(MilliMeter(75.00))
    @test GearModule(pd, nt) ≈ mo
  end


end