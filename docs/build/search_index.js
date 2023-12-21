var documenterSearchIndex = {"docs":
[{"location":"#UnitTypes.jl","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"This package provides physical units as Julia types.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"julia> using UnitTypes\n\njulia> x = Meter(3)\n3m\n\njulia> typeof(x)\nMeter\n\njulia> typeof(x) <: AbstractLength\ntrue\n\njulia> typeof(x) <: AbstractCapacitance\nfalse","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"This allows you to easily write functions with arguments restricted to variables having certain types.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"julia> function goFaster(a::T) where T<:AbstractAcceleration end","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"This leads to correctness and very clear error messages.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"julia> v = MeterPerSecond(3)\n3m/s\n\njulia> goFaster(v)\nERROR: MethodError: no method matching goFaster(::MeterPerSecond)\n\nClosest candidates are:\n  goFaster(::T) where T<:AbstractAcceleration","category":"page"},{"location":"#Introducing-new-types","page":"UnitTypes.jl","title":"Introducing new types","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Macros are used to introduce and create relationships around new types:","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"@makeBaseMeasure Torque NewtonMeter \"N*m\" - introduces a new basic Measure like Meter for Length or Meter3 Volume,\n@deriveMeasure NewtonMeter(1) = MilliNewtonMeter(1000) \"mN*m - introduces a new name for a Measure, often a prefix like Millimeter or an alternate name like Inch, \n@makeDimension Diameter Meter - creates a Dimension, which is a Measure in some particular context, as diameter, radius, and circumference all refer to lengths of a circle.","category":"page"},{"location":"#Design","page":"UnitTypes.jl","title":"Design","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"UnitTypes introduces an abstract type hierarchy of:","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"AbstractMeasure","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Meter, Millimeter, ..., MeterPerSecond, MeterPerSecond2, ... See src/SIDerived.jl Inch, Foot, Mile, ..., See src/Imperial.jl","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"AbstractDimension - src/Dimension.jl","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"AbstractDiameter, AbstractRadius, ... AbstractDuration, ...","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"See src/typeTree.txt for a full list of the pre-defined types.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Please open an issue or PR to add more units.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"As said, the idea is that a Measure is some quantity bearing units, while a Dimension is some context-specific application of a Measure. Within a Dimension multiple Measures may logically be used as long as they are dimensionally consistent. For instance, a circle may be described by its radius, diameter, or circumference, concepts that can be interchangeably converted, using any Measure of extent (<:AbstractLength). A function creating a circle can then internally store radii while accepting Radius, Diameter, or Circumference arguments, as the type system provides conversion between the argument and the function's internal convention.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Concrete Dimensions look like","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"struct Diameter{T <: AbstractLength } <: AbstractDimension\n  value::T\nend","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"and a concrete Measure is represented by","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"struct Meter <: AbstractLength\n  value::Number\n  toBase::Number\n  unit::String\nend","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Measures.jl and the macros define the necessary convert()s and other operators. Please open an issue with a minimal working example if you run into conversion errors.","category":"page"},{"location":"#Logical-operations","page":"UnitTypes.jl","title":"Logical operations","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Using units correctly requires distinguishing between valid and invalid operations, which in some cases means not allowing convenient operations. Inches can be added, as can inch and millimeter, but only when computing area does inch*inch make sense. Inch * 3 is convenient while 3 / Inch is unlikely to be desirable.","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"With use and issues, these coherence rules will become more clear and explained by example.","category":"page"},{"location":"#Comparison-with-other-packages","page":"UnitTypes.jl","title":"Comparison with other packages","text":"","category":"section"},{"location":"#Unitful.jl","page":"UnitTypes.jl","title":"Unitful.jl","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Unitful leverages parametric types to store units, giving flexibility at the cost of compile-time type uncertainty. It's two major limitations are the avoidance of angular measures, as they are not first-class but rather ratios, and rather lengthy type unions that clutter outputs, especially on error:","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"julia> function goSlower(x::T) where T<:Unitful.Acceleration end\ngoSlower (generic function with 1 method)\n\njulia> a = 1u\"mm\"\n\njulia> goSlower(a)\nERROR: MethodError: no method matching goSlower(::Quantity{Int64, 𝐋 , Unitful.FreeUnits{(mm,), 𝐋 , nothing}})\n\nClosest candidates are:\n  goSlower(::T) where T<:(Union{Quantity{T, 𝐋 𝐓^-2, U}, Level{L, S, Quantity{T, 𝐋 𝐓^-2, U}} where {L, S}} where {T, U}) ","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"As Unitful is the dominant unit package and has wide use and support, we provide a separate package ExchangeUnitful to enable interoperation with Unitful.","category":"page"},{"location":"#DynamicQuantities.jl","page":"UnitTypes.jl","title":"DynamicQuantities.jl","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"DynamicQuantities is newer and faster than Unitful because it \"defines a simple statically-typed Quantity type for storing physical units.\" It does this by storing the exponents on the basic units, allowing any unit traceable to SI to be used. But this performant representation hurts readability, and while the unit representation may be able to be hidden behind overrides of show(), Julia is designed for types to be read and manipulated directly by users.","category":"page"},{"location":"#UnitTypes.jl-2","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"In the presence of Julia's type-heavy UI, these two, good attempts feel misdirected and motivate this package's literal typing of units. The limitation is that UnitTypes does not have a catch-all unit representation. Only units that have been defined by one of the macros may be represented, and complex units may need to have additional methods written to correctly convert between units, ie Celsius to Fahrenheit. See SIDerived.jl and Imperial.jl for examples.","category":"page"},{"location":"#Docs","page":"UnitTypes.jl","title":"Docs","text":"","category":"section"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"CurrentModule=UnitTypes","category":"page"},{"location":"","page":"UnitTypes.jl","title":"UnitTypes.jl","text":"Modules=[UnitTypes]","category":"page"},{"location":"#UnitTypes.APerM","page":"UnitTypes.jl","title":"UnitTypes.APerM","text":"This UnitType represents units of APerM.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.APerM2","page":"UnitTypes.jl","title":"UnitTypes.APerM2","text":"This UnitType represents units of APerM2.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Acre","page":"UnitTypes.jl","title":"UnitTypes.Acre","text":"UnitType Acre is derived from Meter2, related by 4046.873, with supertype AbstractArea, and symbol [ac].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Ampere","page":"UnitTypes.jl","title":"UnitTypes.Ampere","text":"This UnitType represents units of Ampere.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Candela","page":"UnitTypes.jl","title":"UnitTypes.Candela","text":"This UnitType represents units of Candela.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.CentiMeter","page":"UnitTypes.jl","title":"UnitTypes.CentiMeter","text":"UnitType CentiMeter is derived from Meter, related by 0.01, with supertype AbstractLength, and symbol [cm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Coulomb","page":"UnitTypes.jl","title":"UnitTypes.Coulomb","text":"This UnitType represents units of Coulomb.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Degree","page":"UnitTypes.jl","title":"UnitTypes.Degree","text":"UnitType Degree is derived from Radian, related by 0.017453292519943295, with supertype AbstractAngle, and symbol [°].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Farad","page":"UnitTypes.jl","title":"UnitTypes.Farad","text":"This UnitType represents units of Farad.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Femtometer","page":"UnitTypes.jl","title":"UnitTypes.Femtometer","text":"UnitType Femtometer is derived from Meter, related by 1.0e-15, with supertype AbstractLength, and symbol [fm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.FluidOunce","page":"UnitTypes.jl","title":"UnitTypes.FluidOunce","text":"UnitType FluidOunce is derived from Meter3, related by 0.0284130625, with supertype AbstractVolume, and symbol [floz].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Foot","page":"UnitTypes.jl","title":"UnitTypes.Foot","text":"UnitType Foot is derived from Meter, related by 0.30479999999999996, with supertype AbstractLength, and symbol [ft].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.FootPerSecond","page":"UnitTypes.jl","title":"UnitTypes.FootPerSecond","text":"UnitType FootPerSecond is derived from MeterPerSecond, related by 0.30479999999999996, with supertype AbstractVelocity, and symbol [ft/s].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Gallon","page":"UnitTypes.jl","title":"UnitTypes.Gallon","text":"UnitType Gallon is derived from Meter3, related by 4.54609, with supertype AbstractVolume, and symbol [gal].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Henry","page":"UnitTypes.jl","title":"UnitTypes.Henry","text":"This UnitType represents units of Henry.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Hertz","page":"UnitTypes.jl","title":"UnitTypes.Hertz","text":"This UnitType represents units of Hertz.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Inch","page":"UnitTypes.jl","title":"UnitTypes.Inch","text":"UnitType Inch is derived from Meter, related by 0.0254, with supertype AbstractLength, and symbol [in].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Kelvin","page":"UnitTypes.jl","title":"UnitTypes.Kelvin","text":"This UnitType represents units of Kelvin.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KgPerM2","page":"UnitTypes.jl","title":"UnitTypes.KgPerM2","text":"This UnitType represents units of KgPerM2.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KgPerM3","page":"UnitTypes.jl","title":"UnitTypes.KgPerM3","text":"This UnitType represents units of KgPerM3.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KiloGram","page":"UnitTypes.jl","title":"UnitTypes.KiloGram","text":"This UnitType represents units of KiloGram.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KiloMeter","page":"UnitTypes.jl","title":"UnitTypes.KiloMeter","text":"UnitType KiloMeter is derived from Meter, related by 1000.0, with supertype AbstractLength, and symbol [km].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KiloNewton","page":"UnitTypes.jl","title":"UnitTypes.KiloNewton","text":"UnitType KiloNewton is derived from Newton, related by 1000.0, with supertype AbstractForce, and symbol [kN].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KiloOhm","page":"UnitTypes.jl","title":"UnitTypes.KiloOhm","text":"UnitType KiloOhm is derived from Ohm, related by 1000.0, with supertype AbstractResistance, and symbol [kΩ].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.KiloVolt","page":"UnitTypes.jl","title":"UnitTypes.KiloVolt","text":"UnitType KiloVolt is derived from Volt, related by 1000.0, with supertype AbstractElectricPotential, and symbol [KV].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Liter","page":"UnitTypes.jl","title":"UnitTypes.Liter","text":"UnitType Liter is derived from Meter3, related by 0.001, with supertype AbstractVolume, and symbol [L].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Lumen","page":"UnitTypes.jl","title":"UnitTypes.Lumen","text":"This UnitType represents units of Lumen.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Lux","page":"UnitTypes.jl","title":"UnitTypes.Lux","text":"This UnitType represents units of Lux.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.M3PerKg","page":"UnitTypes.jl","title":"UnitTypes.M3PerKg","text":"This UnitType represents units of M3PerKg.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MegaOhm","page":"UnitTypes.jl","title":"UnitTypes.MegaOhm","text":"UnitType MegaOhm is derived from Ohm, related by 1.0e6, with supertype AbstractResistance, and symbol [MΩ].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Meter","page":"UnitTypes.jl","title":"UnitTypes.Meter","text":"This UnitType represents units of Meter.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Meter2","page":"UnitTypes.jl","title":"UnitTypes.Meter2","text":"This UnitType represents units of Meter2.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Meter3","page":"UnitTypes.jl","title":"UnitTypes.Meter3","text":"This UnitType represents units of Meter3.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MeterPerSecond","page":"UnitTypes.jl","title":"UnitTypes.MeterPerSecond","text":"This UnitType represents units of MeterPerSecond.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MeterPerSecond2","page":"UnitTypes.jl","title":"UnitTypes.MeterPerSecond2","text":"This UnitType represents units of MeterPerSecond2.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MicroFarad","page":"UnitTypes.jl","title":"UnitTypes.MicroFarad","text":"UnitType MicroFarad is derived from Farad, related by 1.0e-6, with supertype AbstractCapacitance, and symbol [μF].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MicroMeter","page":"UnitTypes.jl","title":"UnitTypes.MicroMeter","text":"UnitType MicroMeter is derived from Meter, related by 1.0e-6, with supertype AbstractLength, and symbol [μm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Mile","page":"UnitTypes.jl","title":"UnitTypes.Mile","text":"UnitType Mile is derived from Meter, related by 1609.3439999999998, with supertype AbstractLength, and symbol [mi].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliFarad","page":"UnitTypes.jl","title":"UnitTypes.MilliFarad","text":"UnitType MilliFarad is derived from Farad, related by 0.001, with supertype AbstractCapacitance, and symbol [mF].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliHenry","page":"UnitTypes.jl","title":"UnitTypes.MilliHenry","text":"UnitType MilliHenry is derived from Henry, related by 0.001, with supertype AbstractInductance, and symbol [mH].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliLiter","page":"UnitTypes.jl","title":"UnitTypes.MilliLiter","text":"UnitType MilliLiter is derived from Liter, related by 0.001, with supertype AbstractVolume, and symbol [mL].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliMeter","page":"UnitTypes.jl","title":"UnitTypes.MilliMeter","text":"UnitType MilliMeter is derived from Meter, related by 0.001, with supertype AbstractLength, and symbol [mm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliNewton","page":"UnitTypes.jl","title":"UnitTypes.MilliNewton","text":"UnitType MilliNewton is derived from Newton, related by 0.001, with supertype AbstractForce, and symbol [mN].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliNewtonMeter","page":"UnitTypes.jl","title":"UnitTypes.MilliNewtonMeter","text":"UnitType MilliNewtonMeter is derived from NewtonMeter, related by 0.001, with supertype AbstractTorque, and symbol [mN*m].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.MilliOhm","page":"UnitTypes.jl","title":"UnitTypes.MilliOhm","text":"UnitType MilliOhm is derived from Ohm, related by 0.001, with supertype AbstractResistance, and symbol [Ω].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.NanoFarad","page":"UnitTypes.jl","title":"UnitTypes.NanoFarad","text":"UnitType NanoFarad is derived from Farad, related by 1.0e-9, with supertype AbstractCapacitance, and symbol [nF].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.NanoMeter","page":"UnitTypes.jl","title":"UnitTypes.NanoMeter","text":"UnitType NanoMeter is derived from Meter, related by 1.0e-9, with supertype AbstractLength, and symbol [nm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.NauticalMile","page":"UnitTypes.jl","title":"UnitTypes.NauticalMile","text":"UnitType NauticalMile is derived from Meter, related by 1852.0, with supertype AbstractLength, and symbol [nmi].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Newton","page":"UnitTypes.jl","title":"UnitTypes.Newton","text":"This UnitType represents units of Newton.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.NewtonMeter","page":"UnitTypes.jl","title":"UnitTypes.NewtonMeter","text":"This UnitType represents units of NewtonMeter.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.NewtonMilliMeter","page":"UnitTypes.jl","title":"UnitTypes.NewtonMilliMeter","text":"UnitType NewtonMilliMeter is derived from NewtonMeter, related by 0.001, with supertype AbstractTorque, and symbol [N*mm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Ohm","page":"UnitTypes.jl","title":"UnitTypes.Ohm","text":"This UnitType represents units of Ohm.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Ounce","page":"UnitTypes.jl","title":"UnitTypes.Ounce","text":"UnitType Ounce is derived from KiloGram, related by 0.028349523125, with supertype AbstractMass, and symbol [oz].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Pascal","page":"UnitTypes.jl","title":"UnitTypes.Pascal","text":"This UnitType represents units of Pascal.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.PerSecond","page":"UnitTypes.jl","title":"UnitTypes.PerSecond","text":"UnitType PerSecond is derived from Hertz, related by 1.0, with supertype AbstractFrequency, and symbol [s^-1].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.PicoFarad","page":"UnitTypes.jl","title":"UnitTypes.PicoFarad","text":"UnitType PicoFarad is derived from Farad, related by 1.0e-12, with supertype AbstractCapacitance, and symbol [pF].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.PicoMeter","page":"UnitTypes.jl","title":"UnitTypes.PicoMeter","text":"UnitType PicoMeter is derived from Meter, related by 1.0e-12, with supertype AbstractLength, and symbol [pm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Pint","page":"UnitTypes.jl","title":"UnitTypes.Pint","text":"UnitType Pint is derived from Meter3, related by 0.56826126, with supertype AbstractVolume, and symbol [pt].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.PoundMass","page":"UnitTypes.jl","title":"UnitTypes.PoundMass","text":"UnitType PoundMass is derived from KiloGram, related by 0.45359237, with supertype AbstractMass, and symbol [lbm].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Quart","page":"UnitTypes.jl","title":"UnitTypes.Quart","text":"UnitType Quart is derived from Meter3, related by 1.1365225, with supertype AbstractVolume, and symbol [qt].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Radian","page":"UnitTypes.jl","title":"UnitTypes.Radian","text":"This UnitType represents units of Radian.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Second","page":"UnitTypes.jl","title":"UnitTypes.Second","text":"This UnitType represents units of Second.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Siemens","page":"UnitTypes.jl","title":"UnitTypes.Siemens","text":"This UnitType represents units of Siemens.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Slug","page":"UnitTypes.jl","title":"UnitTypes.Slug","text":"UnitType Slug is derived from KiloGram, related by 14.59390294, with supertype AbstractMass, and symbol [slug].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.SquareFoot","page":"UnitTypes.jl","title":"UnitTypes.SquareFoot","text":"UnitType SquareFoot is derived from Meter2, related by 0.09290303999999998, with supertype AbstractArea, and symbol [sqft].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.SquareMile","page":"UnitTypes.jl","title":"UnitTypes.SquareMile","text":"UnitType SquareMile is derived from Meter2, related by 2.5899881103359996e6, with supertype AbstractArea, and symbol [sqmi].\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Tesla","page":"UnitTypes.jl","title":"UnitTypes.Tesla","text":"This UnitType represents units of Tesla.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Volt","page":"UnitTypes.jl","title":"UnitTypes.Volt","text":"This UnitType represents units of Volt.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Watt","page":"UnitTypes.jl","title":"UnitTypes.Watt","text":"This UnitType represents units of Watt.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Weber","page":"UnitTypes.jl","title":"UnitTypes.Weber","text":"This UnitType represents units of Weber.\n\n\n\n\n\n","category":"type"},{"location":"#UnitTypes.Yard","page":"UnitTypes.jl","title":"UnitTypes.Yard","text":"UnitType Yard is derived from Meter, related by 0.9143999999999999, with supertype AbstractLength, and symbol [yd].\n\n\n\n\n\n","category":"type"},{"location":"#Base.show-Union{Tuple{T}, Tuple{IO, T}} where T<:AbstractDimension","page":"UnitTypes.jl","title":"Base.show","text":"@show functionality for Dimensions via dimension2String().\n\n\n\n\n\n","category":"method"},{"location":"#UnitTypes.dimension2String-Tuple{T} where T<:AbstractDimension","page":"UnitTypes.jl","title":"UnitTypes.dimension2String","text":"dimension2String(c::AbstractDimension) -> String\n\n\nReturns a string representing dimenson c with format Module.DimensionName(value unit).\n\n\n\n\n\n","category":"method"},{"location":"#UnitTypes.@makeDimension-Tuple{Any, Any}","page":"UnitTypes.jl","title":"UnitTypes.@makeDimension","text":"Make a new dimension dimName of measure; also creates 'AbstractdimName'\n\n```julia     @makeDimension Diameter Meter \n\nd = Diameter(MilliMeter(3.4))\nr = Radius(d)\n\n```\n\n\n\n\n\n","category":"macro"},{"location":"#UnitTypes.@relateDimensions-Tuple{Any}","page":"UnitTypes.jl","title":"UnitTypes.@relateDimensions","text":"Defines various Base. functions that facilitate the given relationship.   All types must already be defined and written in the form type1 = factor * type2, as in:   julia     @relateDimensions Diameter = 2.0*Radius\n\n\n\n\n\n","category":"macro"}]
}