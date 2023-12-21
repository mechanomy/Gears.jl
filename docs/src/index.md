# Gears.jl
This package provides basic functions for working with and plotting spur gears.


* functions for calculating involute tooth profiles
* working with gear parameters and handbook equations
* gear standards
* plotting and other geometric outputs

As these encompass a lot of terrain, the aim is not to be comprehensive but immediately useful, supporting quick prototyping of gears especially through 3D printing.

[Repo](https://github.com/mechanomy/Gears.jl)

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

## Functions

```@meta
CurrentModule=UnitTypes
```

```@autodocs
Modules=[UnitTypes]
```

