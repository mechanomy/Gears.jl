module Gears
  # References:
  # 2004_oberg ~p1599
  # 1989_Hindhede_SpurGearFundamentals - has good diagrams
  # 2012_Dooner_KinematicGeometryOfGearing ch2 - relates tooth width to involute start angle

  using UnitTypes
  using TestItems
  using DocStringExtensions
  using Reexport

  abstract type AbstractGear end
  export AbstractGear

  include("GearDimensions.jl")
  @reexport using .GearDimensions

  include("GearANSI.jl")
  include("InvoluteTooth.jl")

end #Gears

