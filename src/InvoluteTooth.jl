module InvoluteTooth # this submodule provides functions for drawing gear involutes
  using TestItems
  using DocStringExtensions
  using GLMakie

  using Roots
  using LaTeXStrings
  using Printf
  using UnitTypes
  using ..Gears

  # conventions:
  # An involute begins at a point on the base circle at circular angle alpha.
  # It is parameterized by angle theta (>alpha) above the base radius.
  # Angle theta is geometrically interpreted as the point on the base circle whose tangent line is also perpendicular to the involute, as shown in plotInvoluteConstruction.

  """
    Calculates the x coordinate of a point along the involute
    Note that the handedness of the involute is determined by th-al: if th-al>0 a right hand involute, th-al<0 is a left hand
  """
  ix(bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle) = (bd.measure/2*( cos(th) + convert(Radian,th-al).value*sin(th) ) )# returns an AbstractLength...

  """
    Calculates the y coordinate of a point along the involute
  """
  iy(bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle) = (bd.measure/2*( sin(th) - convert(Radian,th-al).value*cos(th) ))
  @testitem "involute calcs" begin
    using UnitTypes
    @test Gears.InvoluteTooth.ix(BaseDiameter(Meter(2)), Radian(1), Radian(2)) ≈ Meter( 1*(cos(2) + (2-1)*sin(2)) ) #check type correctness
    @test Gears.InvoluteTooth.iy(BaseDiameter(Meter(2)), Radian(1), Radian(2)) ≈ Meter( 1*(sin(2) - (2-1)*cos(2)) )
  end

  """
    Gamma is the angle of the line of symmetry of `iTooth` for a complete gear with `nTeeth`
  """
  gamma(nTeeth::Int, iTooth::Int) = Radian(2*π/nTeeth*iTooth) #gamma is the angle of the line of the tooth center

  """
    Calculates the angular tooth width at the pitch diameter = pi/(2*Nteeth).
    See Dooner#61.
  """
  toothAngularWidthAtPitchDiameter(g::AbstractGear) = Radian(π/2/g.nTeeth)

  """
    Finds the angle when the involute has risen `toothHeight` above `r`
    
    $TYPEDSIGNATURES

    Gamma is the angle of the tooth's line of symmetry
    Alpha is the angle of the involute's root on the base circle
    toothHeight is the desired radial height of the tooth beyond the base circle
  """
  function findInvoluteAngleAtRadius(;bd::BaseDiameter, gm::AbstractAngle, al::AbstractAngle, rDesired::AbstractLength)
    ix2(th) = toBaseFloat(ix(bd,al,th))*toBaseFloat(ix(bd,al,th))
    iy2(th) = toBaseFloat(iy(bd,al,th))*toBaseFloat(iy(bd,al,th))
    rinv(th) = sqrt(ix2(th) + iy2(th))
    fal(th) = toBaseFloat(rDesired) - rinv(Radian(th))
    thDesired = Radian(find_zero(fal, convert(Radian,gm).value + convert(Radian,gm-al).value*2))
    return thDesired
  end
  @testitem "findInvoluteAngleAtRadius" begin
    using UnitTypes

    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    alBase = Radian(1)
    gm = Radian(1.2)

    @test isapprox( Gears.InvoluteTooth.findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=g.outside.measure/2), Radian(1.445), atol=1e-3)
  end

  """
    Returns (x,y) positions tracing the gear teeth for gear `g`. 
    (x,y) is a vector of Measures which <:AbstractLength.
  """
  function getToothProfilePoints(g::G where G<:Gears.AbstractGear; nPerTooth::Int=100)
    nPerInvolute = convert(Int64,round(nPerTooth/2)) #this may be better defined by physical spacing...

    xs = Meter.(zeros(g.nTeeth*nPerInvolute*2+1))
    ys = Meter.(zeros(g.nTeeth*nPerInvolute*2+1))
    for i in 1:round(g.nTeeth)
      gm=gamma(g.nTeeth, i) 

      psi = Radian(acos(g.base.measure.value/g.pitch.measure.value))
      vpsi = Radian(tan(psi))-psi

      # tooth frontside
      alPitch = gm - toothAngularWidthAtPitchDiameter(g) #cf fig2.4 2012Dooner
      alBase = alPitch - vpsi 
      thOutside = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=g.outside.measure/2) # this is the angle of the involute tip at the outside diameter
      thRoot = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=g.root.measure/2)
      ths = Radian.(LinRange(convert(Radian,thRoot).value, convert(Radian,thOutside).value, nPerInvolute))
      for j = 1:nPerInvolute
        xs[(i-1)*2*nPerInvolute + j] = Gears.InvoluteTooth.ix(g.base, alBase, ths[j])
        ys[(i-1)*2*nPerInvolute + j] = Gears.InvoluteTooth.iy(g.base, alBase, ths[j])
      end

      # tooth backside
      btPitch = gm + toothAngularWidthAtPitchDiameter(g) #cf fig2.4 2012Dooner
      btBase = btPitch + vpsi 
      thOutside = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=btBase, rDesired=g.outside.measure/2) # this is the angle of the involute tip at the outside diameter
      thRoot = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=btBase, rDesired=g.root.measure/2)
      ths = Radian.(LinRange(convert(Radian,thOutside).value, convert(Radian,thRoot).value, nPerInvolute))
      for j = 1:nPerInvolute
        xs[(i-1)*2*nPerInvolute + nPerInvolute + j] = Gears.InvoluteTooth.ix(g.base, btBase, ths[j])
        ys[(i-1)*2*nPerInvolute + nPerInvolute + j] = Gears.InvoluteTooth.iy(g.base, btBase, ths[j])
      end

    end
    xs[length(xs)] = xs[1] # close the gear
    ys[length(ys)] = ys[1]

    return (xs,ys)
  end
  @testitem "getToothProfilePoints" begin
    using UnitTypes
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    (xs,ys) = Gears.InvoluteTooth.getToothProfilePoints(g, nPerTooth=20)

    # inspect the plots produce by plotInvoluteConstruction and plotGearTeeth, if they look correct hard-code a few values
    # @show xs[1:3] last(xs)
    # @show ys[1:3] last(ys)
    @test isapprox( xs[1], Meter(.03566), atol=1e-4 )
    @test isapprox( xs[2], Meter(.03585), atol=1e-4 ) 
    @test isapprox( xs[3], Meter(.03606), atol=1e-4 ) 
    @test isapprox( last(xs), Meter(.03566), atol=1e-4 ) 
    @test isapprox( ys[1], Meter(.00201), atol=1e-4 )
    @test isapprox( ys[2], Meter(.00207), atol=1e-4 ) 
    @test isapprox( ys[3], Meter(.00214), atol=1e-4 ) 
    @test isapprox( last(ys), Meter(.00201), atol=1e-4 ) 
  end

  """
    Writes a file of points tracing the outer edge of the gear teeth.

    $TYPEDSIGNATURES

    * `g` the Gear to write
    * `fileName` optional file name, if not provided will be named "gearProfilePoints_diametralPitch_nTeeth"*fileExtension
    * `fileExtension` as below
    * `nPerTooth` the number of points to write per tooth
    * `unitType` the unit to save the points in from UnitTypes

    The file format is given by the file extension in `fileName` or `fileExtension`.
    Suported file types are:
    * tsv==txt = tab-separated values
    * csv = comma-separated values
    * sldcrv = Solidworks Curve format having rows of Xin Yin Zin, where XYZ are the float positions and 'in' the linear unit. As the tooth profile is XY, Z = 0.

    ```julia
      using Gears
      using UnitTypes
      g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) )
      Gears.InvoluteTooth.writeToothProfilePoints(g, fileName="gear70.csv", unitType=CentiMeter)
      Gears.InvoluteTooth.writeToothProfilePoints(g, fileName="gear70", fileExtension=".txt", unitType=Inch)
    ```
  """
  function writeToothProfilePoints(g::Gears.AbstractGear; fileName::String="", fileExtension::String=splitext(fileName)[2], nPerTooth::Int=100, unitType::Type{T}=AbstractLength) where T<:AbstractLength
    (xs, ys) = getToothProfilePoints(g, nPerTooth=nPerTooth)

    # determine the file format by the given extension
    if fileName=="" 
      if fileExtension==""
        fileExtension = ".txt"
      end
      fileName = "gearProfilePoints_$(g.diametral)_$(g.nTeeth)"*fileExtension
    end

    if splitext(fileName)[2] == ""
      fileName = fileName * fileExtension
    end

    #now emit the points in the format expected by CAD
    zs = xs *0

    # which unit to emit?
    # ys = convert(typeof(xs[1]), ys)
    if unitType == AbstractLength #default
      xs = map(x -> convert(MilliMeter,x), xs)
      ys = map(x -> convert(MilliMeter,x), ys)
    else
      xs = map(x -> convert(unitType,x), xs)
      ys = map(x -> convert(unitType,x), ys)
    end

    un = xs[1].unit

    open(fileName, "w") do f
      if lowercase(fileExtension) == ".tsv" || lowercase(fileExtension) == ".txt"
        write(f, Gears.gear2String(g)*"\n" )
        write(f,"x[$un]\ty[$un]\n")
        for i in eachindex(xs)
          write(f, @sprintf("%3.3f\t%3.3f\n", xs[i].value, ys[i].value))
        end

      elseif lowercase(fileExtension) == ".csv"
        write(f, Gears.gear2String(g)*"\n" )
        write(f,"x[$un],y[$un]\n")
        for i in eachindex(xs)
          write(f, @sprintf("%3.3f,%3.3f\n", xs[i].value, ys[i].value))
        end

      elseif fileExtension == lowercase(".sldcrv")
        for i in eachindex(xs)
          #solidworks format: 25mmTAB0mmTAB0mmCRLF
          write(f, "$(xs[i].value)$un\t$(ys[i].value)$un\t$(zs[i].value)$un\r\n") # inches, choose unit by?
        end
      else
        throw(ArgumentError("File extension [$fileExtension] is not known to Gears.InvoluteTooth.writeToothProfilePoints()"))
      end
    end #open
  end
  @testitem "writeToothProfilePoints to files" begin
    # to test, check that the file exists and the first line after the header is correct
    using UnitTypes
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070

    tDir = mktempdir()
    println("Test folder is $tDir")

    tPath = joinpath(tDir, "testProfilePoints.txt") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath)
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      l2 = readline(fid)
      l3 = readline(fid)
      # @test l3 == "34.787\t-8.106" #as mm
    end

    tPath = joinpath(tDir, "testProfilePoints.tsv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, unitType=CentiMeter)
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      l2 = readline(fid)
      l3 = readline(fid)
      # @test l3 == "3.479\t-0.811" #as cm
    end

    tPath = joinpath(tDir, "testProfilePoints.csv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, fileExtension=".csv", unitType=Inch)
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      l2 = readline(fid)
      l3 = readline(fid)
      # @test l3 == "1.370,-0.319" # as inch
    end

    tPath = joinpath(tDir, "testProfilePoints.tsv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, fileExtension=".sldcrv")
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      # @test l1 == "34.78722899491671mm\t-8.105983123125842mm\t0.0mm"
    end   
  end

  function arcXY(; r::AbstractLength, a0::AbstractAngle, a1::AbstractAngle, n=100)
    ths = LinRange(convert(Radian,a0).value, convert(Radian,a1).value, n) 
    xs = r.value .*cos.(ths)
    ys = r.value .*sin.(ths)
    return (xs,ys)
  end

  """
    Plots the arc between angles
  """
  function plotArc(; axs, r::AbstractLength, a0::AbstractAngle, a1::AbstractAngle, n=100, label="", linecolor=:gray, linestyle=:solid, fontsize=10)
    (xs,ys) = arcXY(r=r, a0=a0, a1=a1, n=n)
    lines!(axs, xs, ys, color=linecolor, linestyle=linestyle)
    if label != ""
      text!(axs, r.value*cos(a0), r.value*sin(a0), text=label, color=linecolor, fontsize=fontsize, align=(:left,:top) ) #place near a0
    end
    return (xs,ys)
  end

  """
    Plots an involute of the circle
  """
  function plotInvolute(; bd::BaseDiameter, al::AbstractAngle, thMax::AbstractAngle, axs=nothing, linecolor=:blue, linestyle=:solid)
    if isnothing(axs)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    ths = Radian.(LinRange(convert(Radian,al).value, convert(Radian,thMax).value, 100))
    xs = zeros(length(ths))
    ys = zeros(length(ths))
    for i = eachindex(ths)
      xs[i] = Gears.InvoluteTooth.ix(bd, al, ths[i]).value #again just take the raw value, waiting for unit-aware graphs
      ys[i] = Gears.InvoluteTooth.iy(bd, al, ths[i]).value
    end

    l = lines!(ax,xs,ys, color=linecolor, linestyle=linestyle)

    if isnothing(axs)
      display(GLMakie.Screen(), fig)
    end

    return l
  end

  """
    Plots two lines that explain the construction of the involute at the given angle.
  """
  function plotInvoluteConstruction(; bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle, axs=nothing, linecolor=:gray, linestyle=:solid)
    if isnothing(axs)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end
    r = bd.measure.value/2 # take in the units given for now; on rewrite I should establish a unit for the graph and convert all

    xs = [0, r*cos(th), ix(bd,al,th).value]
    ys = [0, r*sin(th), iy(bd,al,th).value]
    l = lines!(ax, xs,ys, color=linecolor, linestyle=linestyle)

    if isnothing(axs)
      display(GLMakie.Screen(), fig)
    end

    return l
  end

  function plotAngleArc(; axs, r::AbstractLength, aMax::AbstractAngle, aMin::AbstractAngle=Radian(0), label="angleArc@$aMax", aLabel::AbstractAngle=aMax/2, linecolor=:gray, linestyle=:solid, fontsize=10)
    ths = Radian.(LinRange(convert(Radian,aMin).value, convert(Radian,aMax).value, 100))
    xs = r.value .* cos.(ths)
    ys = r.value .* sin.(ths)
    s = lines!(axs, xs, ys, color=linecolor, linestyle=linestyle)
    text!(axs, r.value*cos(aLabel), r.value*sin(aLabel), text=label, color=linecolor, fontsize=fontsize, align=(:left,:bottom) )
    return s
  end

  function plotAngleLine(;axs, r::AbstractLength, a::AbstractAngle, rAngle::AbstractLength=r*.8, label="", linecolor=:gray, linestyle=:solid, fontsize=10)
    plotAngleArc(axs=axs, r=rAngle, aMax=a, aMin=Radian(0), label=label )
    lines!(axs, [0,r.value*cos(a)], [0,r.value*sin(a)], color=linecolor, linestyle=linestyle, label=label)
    # if label != ""
    #   text!(axs, r.value/2*cos(a), r.value/2*sin(a), text=label, color=linecolor, fontsize=fontsize, align=(:left,:top) )
    # end
  end

  """
    This function draws an annottated plot showing the key geometric relations for the given Gear `g`.
  """
  function plotInvoluteConstruction(g::AbstractGear, toothNumber=Int(round(g.nTeeth/8)), fig=nothing )
    un = Inch #unit for x and y

    # println("\n>>Gear base[$(g.base)] root[$(g.root)] pitch[$(g.pitch)] outside[$(g.outside)]")
    if isnothing(fig)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      axs = Axis(fig[1,1], xlabel="X[$(un(1).unit)]", ylabel="Y[$(un(1).unit)]", title=string(g), aspect=DataAspect())
    else
      # axs = Axis(fig) #assume that fig has only one axis?
      axs = fig[1,1]
    end
    
    gm = gamma(g.nTeeth, toothNumber) #gamma is the tooth centerline
    gmExtra = 0.2
    plotArc(axs=axs, r=g.base.measure/2, a0=gm-gmExtra, a1=gm+gmExtra, label="base", linecolor=:black)
    plotArc(axs=axs, r=g.root.measure/2, a0=gm-gmExtra, a1=gm+gmExtra, label="root", linecolor=:blue)
    plotArc(axs=axs, r=g.pitch.measure/2, a0=gm-gmExtra, a1=gm+gmExtra, label="pitch", linecolor=:green)
    plotArc(axs=axs, r=g.outside.measure/2, a0=gm-gmExtra, a1=gm+gmExtra, label="outside", linecolor=:red) # = addendum circle

    ###############################################################################################
    # the first involute on the <gm side of the tooth, located on the base circle by alpha
    # the alpha angles locate the point on the base circle where a tangent line is perpendicular to the involute
    psi = Radian(acos(toBaseFloat(g.base.measure)/toBaseFloat(g.pitch.measure)))
    vpsi = Radian(tan(psi))-psi

    alPitch = gm - toothAngularWidthAtPitchDiameter(g) #cf fig2.4 2012Dooner
    alBase = alPitch - vpsi 
    thPitch = alPitch + psi 
    thOutside = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=alBase, rDesired=g.outside.measure/2) # this is the angle of the involute tip at the outside diameter

    plotInvolute(axs=axs, bd=g.base, al=alBase, thMax=thOutside, linecolor=:orange)

    ro = g.outside.measure/2 #a reference for spacing the labels
    plotAngleLine(axs=axs, r=g.outside.measure/2, rAngle=ro*.60, a=gm, linecolor=:red, label="γ")
    plotAngleLine(axs=axs, r=g.pitch.measure/2, rAngle=ro*.57, a=alPitch, linecolor=:gray50, label="αPitch" )
    plotAngleLine(axs=axs, r=g.base.measure/2, rAngle=ro*.54, a=alBase, linecolor=:gray50, label="αBase" )
    plotArc(axs=axs, r=ro*.56, a0=alPitch, a1=thPitch, label="ψ", linecolor=:aqua)
    plotArc(axs=axs, r=ro*.55, a0=alBase, a1=alPitch, label="vψ", linecolor=:lightblue)

    plotInvoluteConstruction(axs=axs, bd=g.base, al=alBase, th=thPitch, linecolor=:orange)

    ###############################################################################################
    #now plot the other side of the tooth
    btPitch = gm + toothAngularWidthAtPitchDiameter(g)
    btBase = btPitch + vpsi
    thOutside = findInvoluteAngleAtRadius(bd=g.base, gm=gm, al=btBase, rDesired=g.outside.measure/2) # this is the angle of the involute tip at the outside diameter

    plotAngleLine(axs=axs, r=g.pitch.measure/2, rAngle=ro*.48, a=btPitch, linecolor=:gray80, label="βPitch" )
    plotAngleLine(axs=axs, r=g.base.measure/2, rAngle=ro*.51, a=btBase, linecolor=:gray80, label="βBase" )

    plotInvolute(axs=axs, bd=g.base, al=btBase, thMax=thOutside, linecolor=:magenta)
    thBtPitch = btPitch - psi 
    plotInvoluteConstruction(axs=axs, bd=g.base, al=btBase, th=thBtPitch, linecolor=:magenta)

    rt = ro*1.1
    (xs,ys) = plotArc(axs=axs, r=rt, a0=alPitch, a1=btPitch, label="t = π/N", linecolor=:gray70)
    lines!(axs, [xs[1],ix(g.base,alBase,thPitch).value], [ys[1],iy(g.base,alBase,thPitch).value], color=:gray70)
    lines!(axs, [last(xs),ix(g.base,btBase,thBtPitch).value], [last(ys),iy(g.base,btBase,thBtPitch).value], color=:gray70)

    return fig
  end
  @testitem "plotInvoluteConstruction" begin
    # using UnitTypes
    # using GLMakie
    # GLMakie.closeall()
    # GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config
    # g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    # Gears.InvoluteTooth.plotInvoluteConstruction(g)

    # save("w:/sync/mechgits/julia/Gears/involuteConstruction.png", fig) # GLMakie can only save png
    # DataInspector(fig, transparency=true, backgroundcolor=RGBAf(1,1,1,0.5)) # https://docs.makie.org/stable/explanations/inspector/index.html
    # display(GLMakie.Screen(), fig) # note the window only lasts as long as the julia session
    @test true
  end

  """
    Plots both sides to form a tooth
    r is the root diameter
    gamma is the center angle of the tooth
    delta is the angular width of the tooth at the root diameter
    height is the tooth height, or the difference between the outside and root diameters
  """
  function plotGearTeeth( g::AbstractGear, fig=nothing)
    un = Inch #unit for x and y

    if isnothing(fig)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      axs = Axis(fig[1,1], xlabel="X[$(un(1).unit)]", ylabel="Y[$(un(1).unit)]", title=string(g), aspect=DataAspect())
    else
      # axs = Axis(fig) #assume that fig has only one axis?
      axs = fig[1,1]
    end

    (xs,ys) = getToothProfilePoints(g)

    xs = map(x -> convert(un,x).value, xs)
    ys = map(x -> convert(un,x).value, ys)
    lines!(axs, xs, ys)

    return fig
  end 
  @testitem "plotGearTeeth" begin
    # using UnitTypes
    # using GLMakie
    # GLMakie.closeall()
    # GLMakie.activate!(decorated=true, focus_on_show=true) # https://docs.makie.org/stable/explanations/backends/glmakie/index.html#activation_and_screen_config
    # g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070

    # fig = Gears.InvoluteTooth.plotGearTeeth(g)
    # # fig = Gears.InvoluteTooth.plotInvoluteConstruction(g, 10, fig)
    # # save("w:/sync/mechgits/julia/Gears/gear.png", fig) # GLMakie can only save png

    # DataInspector(fig, transparency=true, backgroundcolor=RGBAf(1,1,1,0.5)) # https://docs.makie.org/stable/explanations/inspector/index.html
    # display(GLMakie.Screen(), fig) # note the window only lasts as long as the julia session
    @test true
  end

end #InvoluteTooth
