# Gears.jl
This package models basic spur gears, including:
* functions for calculating involute tooth profiles
* working with gear parameters and handbook equations
* gear standards
* plotting and other geometric outputs

As these encompass a lot of terrain, the aim is not to be comprehensive but immediately useful, supporting quick prototyping of gears especially through 3D printing.

![](./gear.png)

## Examples

### Plot a gear profile
```julia
using UnitTypes
using GLMakie
g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) 
fig = Gears.InvoluteTooth.plotInvoluteConstruction(g)
fig = Gears.InvoluteTooth.plotGearTeeth(g, fig)
display(fig)
```

### Write gear profile points to a file
```julia
using UnitTypes
g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) 
Gears.InvoluteTooth.writeToothProfilePoints(g, fileName= "testProfilePoints.txt") 
```

## References
* 2012 Dooner "Kinematic Geometry of Gearing"

<!-- Gear Technology has a great archive:

* [Design of Involute Gear Teeth](https://www.geartechnology.com/articles/20828-design-of-involute-gear-teeth)
_1984_Fellows_DesignOfInvoluteGearTeeth.pdf_

* [Determination of gear ratios](https://www.geartechnology.com/articles/20818-determination-of-gear-ratios)
_1984_Orthwein_DeterminationOfGearRatios.pdf_

* [Functions of Gearing and Application of the Involute to Gear Teeth](https://www.geartechnology.com/articles/20819-functions-of-gearing-and-application-of-the-involute-to-gear-teeth)
1984

* [The Design and Manufacture of Machined Plastic Gears](https://www.geartechnology.com/articles/20191-the-design-and-manufacture-of-machined-plastic-gears)
_1985_ChenJuarbe_

* 1985 https://www.geartechnology.com/articles/20103-finding-gear-teeth-ratios

* 1985 https://www.geartechnology.com/articles/20181-gear-inspection-and-chart-interpretation

* [Involumetry](https://www.geartechnology.com/articles/20938-involutometry)
_1988_vanGerpenReece_Involumetry.pdf_

* [Involumetry](https://www.geartechnology.com/articles/20951-involutometry-illustrations)
_1988_vanGerpenReece_InvolumetryIllustrations.pdf_

* [Spur Gear Fundamentals](https://www.geartechnology.com/articles/20954-spur-gear-fundamentals)
_1989_Hindhede_SpurGearFundamentals.pdf_

* [Dudley's Handbook of Practical Gear Design and Manufacture](https://www.google.com/books/edition/Dudley_s_Handbook_of_Practical_Gear_Desi/mLUclpsfTeQC?hl=en)
_2012_Radzevich_DudleysHnadbookOfPracticalGearDesignAndManufacture_ -->


## Copyright
Copyright (c) 2022 Mechanomy LLC

Released under MIT.
