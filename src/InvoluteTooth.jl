module InvoluteTooth # this submodule provides functions for drawing gear involutes
  using TestItems
  using DocStringExtensions
  using GLMakie

  using Roots
  using LaTeXStrings
  using Printf
  using UnitTypes
  using ..Gears

  """
  calculates the x coordinate of a point along the involute
  Note that the handedness of the involute is determined by th-al: if th-al>0 a right hand involute, th-al<0 is a left hand
  """
  ix(bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle) = (Radius(bd.measure/2)*( cos(th) + convert(Radian,th-al).value*sin(th) ) ).measure # returns an AbstractLength...

  """
  calculates the y coordinate of a point along the involute
  """
  iy(bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle) = (Radius(bd.measure/2)*( sin(th) - convert(Radian,th-al).value*cos(th) )).measure
  @testitem "involute calcs" begin
    using UnitTypes
    @test Gears.InvoluteTooth.ix(BaseDiameter(Meter(2)), Radian(1), Radian(2)) ≈ Meter( 1*(cos(2) + (2-1)*sin(2)) ) #check type correctness
    @test Gears.InvoluteTooth.iy(BaseDiameter(Meter(2)), Radian(1), Radian(2)) ≈ Meter( 1*(sin(2) - (2-1)*cos(2)) )
  end


  """
  finds thMax to intersect with the line of symmetry gamma 
  Gamma is the angle of the tooth's line of symmetry
  Alpha is the angle of the involute's root on the base circle
  Theta0 is the guess of the angle corresponding to the intersection of the involute and gamma
  """
  function calcThetaIntersect(; gm::AbstractAngle, al::AbstractAngle, th0::AbstractAngle=gm+(gm-al)*2)
    #calculate thMax to intersect with the gamma line
    fxi(thc,al,gm) = cos(thc)^2 + 2*convert(Radian,thc-al).value*cos(thc)*sin(thc) + (sin(thc)^2-cos(gm)^2)*convert(Radian,thc-al).value^2 - cos(gm)^2
    fyi(thc,al,gm) = sin(thc)^2 - 2*convert(Radian,thc-al).value*cos(thc)*sin(thc) + (cos(thc)^2-sin(gm)^2)*convert(Radian,thc-al).value^2 - sin(gm)^2
    fxz(thc) = fxi(thc, al, gm)^2 + fyi(thc, al, gm)^2 
    return find_zero(fxz, convert(Radian,th0).value, atol=1e-5) 
  end
  @testitem "calcThetaIntersect" begin
    using UnitTypes
    @test isapprox( Gears.InvoluteTooth.calcThetaIntersect(gm=Radian(2), al=Radian(1)), 3.1327, atol=1e-3)
  end

  """
  Finds the angle theta where the involute intersects the root circle, which is where the tooth begins.

  $TYPEDSIGNATURES
    `rd` - the RootDiameter
    `bd` - the BaseDiameter
    `alb` - the angle of the involute at the base diameter
    `al0` - the angle of the involute at the outside diameter
    
  """
  function calcThetaRootIntersection(; rd::RootDiameter, bd::BaseDiameter, alb::AbstractAngle, al0::AbstractAngle )
    fal(thc) = ( rd/2 - bd/2*sqrt(1+convert(Radian,thc-alb).value^2) ).measure.value #just keep the float
    return Radian(find_zero(fal, convert(Radian,al0).value, atol=1e-5))
  end
  @testitem "calcThetaRootIntersection" begin
    using UnitTypes
    g = GearANSI( PitchDiameter(Inch(2.9167)), 70, Degree(20) ) #sdpsi_S1268Z-024A070
    @test isapprox(Gears.InvoluteTooth.calcThetaRootIntersection( rd=g.root, bd=g.base, alb=Radian(1), al0=Radian(0.9) ), 0.7697, atol=1e-3)
  end

  """
  Finds the angle when the involute has risen `toothHeight` above `r`
  
  $TYPEDSIGNATURES

    Gamma is the angle of the tooth's line of symmetry
    Alpha is the angle of the involute's root on the base circle
    toothHeight is the desired radial height of the tooth beyond the base circle
  """
  function calcThetaHeight(; bd::BaseDiameter, gm::AbstractAngle, al::AbstractAngle, toothHeight::AbstractLength)
    htal(th) = toBaseFloat(bd.measure/2) * ( cos(th-gm) + convert(Radian,th-al).value*sin(th-gm) - 1 )
    fal(th) = toBaseFloat(toothHeight) - htal(th)
    thal = Radian(find_zero(fal, convert(Radian,gm).value + convert(Radian,gm-al).value*2))
    return thal
  end
  @testitem "calcThetaHeight" begin
    using UnitTypes
    @test isapprox( Gears.InvoluteTooth.calcThetaHeight(bd=BaseDiameter(Millimeter(192)), gm=Degree(56), al=Radian(0.3), toothHeight=Millimeter(20) ), -0.4870, atol=1e-3)
    @test isapprox( Gears.InvoluteTooth.calcThetaHeight(bd=BaseDiameter(Millimeter(192)), gm=Degree(56), al=Radian(1.27), toothHeight=Millimeter(20) ), 0.4910, atol=1e-3)
    base = BaseDiameter(Millimeter(192))
    outside = OutsideDiameter(Millimeter(200))
    @test isapprox( Gears.InvoluteTooth.calcThetaHeight(bd=base, gm=Degree(56), al=Radian(1.27), toothHeight=(outside-base).measure/2 ), 0.69537, atol=1e-3)
  end

  """
    Returns (x,y) positions tracing the gear teeth for gear `g`. 
    (x,y) is a vector of Measures which <:AbstractLength.

  """
  function getToothProfilePoints(g::G where G<:Gears.AbstractGear; nPerTooth::Int=100)
    nPerInvolute = convert(Int64,round(nPerTooth/2)) #this may be better defined by physical spacing...
    rpd = toBaseFloat(g.pitch.measure)/2
    t = π/(2*g.nTeeth)  # calculating the tooth width t at the pitch diameter = pi/(2*Nteeth) cf Dooner61
  
    xs = Meter.(zeros(g.nTeeth*nPerInvolute*2+1))
    ys = Meter.(zeros(g.nTeeth*nPerInvolute*2+1))
    for i in 1:round(g.nTeeth)
      gm=Radian(2*π/g.nTeeth*(i)) #gamma is the angle of the line of the tooth center

      psi = Radian(acos(g.base.measure.value/g.pitch.measure.value))
      vpsi = Radian(tan(psi))-psi

      # tooth frontside
      alb = gm - Radian(asin(t/2/rpd)) - vpsi
      alo = Gears.InvoluteTooth.calcThetaHeight(bd=g.base, gm=gm, al=alb, toothHeight=(g.outside-g.base).measure/2) # this is the angle of the involute tip at the outside diameter
      thalr = Gears.InvoluteTooth.calcThetaRootIntersection(rd=g.root, bd=g.base, alb=alb, al0=alb)
      ths = Radian.(LinRange(convert(Radian,thalr).value, convert(Radian,alo).value, nPerInvolute))
      for j = 1:nPerInvolute
        xs[(i-1)*2*nPerInvolute + j] = Gears.InvoluteTooth.ix(g.base, alb, ths[j])
        ys[(i-1)*2*nPerInvolute + j] = Gears.InvoluteTooth.iy(g.base, alb, ths[j])
      end

      # tooth backside
      btb = gm + asin(t/2/rpd) + vpsi
      bto = Gears.InvoluteTooth.calcThetaHeight(bd=g.base, gm=gm, al=btb, toothHeight=(g.outside-g.base).measure/2) 
      thbtr = Gears.InvoluteTooth.calcThetaRootIntersection(rd=g.root, bd=g.base, alb=btb, al0=bto)

      ths = Radian.(LinRange(convert(Radian,bto).value, convert(Radian,thbtr).value, nPerInvolute))
      for j = 1:nPerInvolute
        xs[(i-1)*2*nPerInvolute + nPerInvolute + j] = Gears.InvoluteTooth.ix(g.base, btb, ths[j])
        ys[(i-1)*2*nPerInvolute + nPerInvolute + j] = Gears.InvoluteTooth.iy(g.base, btb, ths[j])
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

    # the plot looks correct, so hard-code a few values
    # @show xs[1:3] last(xs)
    # @show ys[1:3] last(ys)
    @test isapprox( xs[1], Meter(0.03479), atol=1e-5 )
    @test isapprox( xs[2], Meter(0.03509), atol=1e-5 ) 
    @test isapprox( xs[3], Meter(0.03543), atol=1e-5 ) 
    @test isapprox( last(xs), Meter(0.03479), atol=1e-5 ) 
    @test isapprox( ys[1], Meter(-0.00811), atol=1e-3 )
    @test isapprox( ys[2], Meter(-0.00810), atol=1e-3 ) 
    @test isapprox( ys[3], Meter(-0.00808), atol=1e-3 ) 
    @test isapprox( last(ys), Meter(-0.00811), atol=1e-3 ) 
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
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName="gear70.csv", unitType=Centimeter)
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
      xs = map(x -> convert(Millimeter,x), xs)
      ys = map(x -> convert(Millimeter,x), ys)
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
      @test l3 == "34.787\t-8.106" #as mm
    end

    tPath = joinpath(tDir, "testProfilePoints.tsv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, unitType=Centimeter)
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      l2 = readline(fid)
      l3 = readline(fid)
      @test l3 == "3.479\t-0.811" #as cm
    end

    tPath = joinpath(tDir, "testProfilePoints.csv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, fileExtension=".csv", unitType=Inch)
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      l2 = readline(fid)
      l3 = readline(fid)
      @test l3 == "1.370,-0.319" # as inch
    end

    tPath = joinpath(tDir, "testProfilePoints.tsv") 
    Gears.InvoluteTooth.writeToothProfilePoints(g, fileName=tPath, fileExtension=".sldcrv")
    @test isfile(tPath)
    open(tPath, "r") do fid
      l1 = readline(fid)
      @test l1 == "34.78722899491671mm\t-8.105983123125842mm\t0.0mm"
    end   
  end


  """
  Plots the base arc under the tooth
  """
  function plotBaseSection(; bd::BaseDiameter, al::AbstractAngle, bt::AbstractAngle, thExtra::AbstractAngle=Radian(0.5), n=100, axs=nothing, linecolor=:gray, linestyle=:solid)
    if isnothing(axs)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    #want to plot only a section of the base circle, that from al to thMax
    ths = LinRange(convert(Radian,al-thExtra).value, convert(Radian,bt+thExtra).value, n) 
    xs = bd.measure.value/2 .*cos.(ths)
    ys = bd.measure.value/2 .*sin.(ths)
    l = lines!(ax, xs, ys, color=linecolor, linestyle=linestyle)

    if isnothing(axs)
      display(GLMakie.Screen(), fig)
    end
    return l
  end

  """
  Plots the line of symmetry gamma
  """
  function plotGamma(; axs, gm::AbstractAngle, bd::BaseDiameter, ht::AbstractLength, linecolor=:green, linestyle=:dash )
    r = bd.measure.value/2 # take in the units given for now; on rewrite I should establish a unit for the graph and convert all
    l = lines!(axs, [0,(r+ht.value)*cos(gm)], [0,(r+ht.value)*sin(gm)], color=linecolor, linestyle=linestyle)
    return l
  end

  """
  Plots the line to the base of the involute
  """
  function plotAlpha(; axs, bd::BaseDiameter, al::AbstractAngle, linecolor=:gray, linestyle=:dash)
    r = bd.measure.value/2 # take in the units given for now; on rewrite I should establish a unit for the graph and convert all
    l = lines!(axs, [0,r*cos(al)], [0,r*sin(al)], color=linecolor, linestyle=linestyle)
    return l
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
  """
  function plotInvoluteConstruction(; bd::BaseDiameter, al::AbstractAngle, th::AbstractAngle, axs=nothing, linecolor=:gray, linestyle=:dash)
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

  """
  """
  function drawAngleArc(; axs, r::AbstractLength, aMax::AbstractAngle, aMin::AbstractAngle=0, label="label", aLabel::AbstractAngle=aMax/2, linecolor=:gray, linestyle=:solid, fontsize=10)
    ths = Radian.(LinRange(convert(Radian,aMin).value, convert(Radian,aMax).value, 100))
    xs = r.value .* cos.(ths)
    ys = r.value .* sin.(ths)

    s = scatterlines!(axs, xs, ys, color=linecolor, linestyle=linestyle, marker=:rtriangle)

    text!(axs, r.value*cos(aLabel), r.value*sin(aLabel), text=label, color=linecolor, fontsize=fontsize, align=(:left,:bottom) )
    return s
  end

  """
  """
  function drawToothTop(;bd::BaseDiameter, al::AbstractAngle, thal::AbstractAngle, bt::AbstractAngle, thbt::AbstractAngle, axs=nothing, linecolor=:orange, linestyle=:dash)
    if isnothing(axs)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect())
    else
      ax = axs
    end

    l = lines!(ax, [ix(bd,al,thal).value, ix(bd,bt,thbt).value], [iy(bd,al,thal).value, iy(bd,bt,thbt).value], color=linecolor, linestyle=linestyle) # just two points = flat, though could assume it's circular at OutsideDiameter..

    if isnothing(axs)
      display(GLMakie.Screen(), fig)
    end
    return l
  end

  """
  Plots both sides to form a tooth
  r is the root diameter
  gamma is the center angle of the tooth
  delta is the angular width of the tooth at the root diameter
  height is the tooth height, or the difference between the outside and root diameters
  """
  function plotTooth(; bd::BaseDiameter, gm::AbstractAngle, dl::AbstractAngle, ht::AbstractLength, axs=nothing)
    # plot construction of both involutes
    al = gm - dl/2
    althMax = calcThetaHeight(bd=bd, gm=gm, al=al, toothHeight=ht)
    bt = gm + dl/2
    btthMax = calcThetaHeight(bd=bd, gm=gm, al=bt, toothHeight=ht)

    if isnothing(axs)
      fig = Figure(backgroundcolor="#bbb", size=(1000,1000))
      ax = Axis(fig[1,1], xlabel="X", ylabel="Y", aspect=DataAspect(), title="Involute tooth with γ=$(Degree(gm).value), δ=$(Degree(dl).value)")
    else
      ax = axs
    end


    lines!(ax, [0,bd.measure.value/2*1.1],[0,0], color=:gray) # draw angle=0 line
    plotBaseSection(axs=ax, bd=bd, al=al, bt=bt, linecolor=:black, linestyle=:solid)
    plotGamma(axs=ax, bd=bd, gm=gm, ht=ht)

    # drawAngleArc(axs=ax, r=40, aMin=0, aMax=gm, aLabel=deg2rad(0), label="γ", color=:gray)
    # savefig("on220826_involutes_baseGamma.svg")
    # drawAngleArc(axs=ax, r=70, aMin=al, aMax=bt, aLabel=deg2rad(35), label=L"\delta", color=:gray)
    # savefig("on220826_involutes_delta.svg")

    plotAlpha(axs=ax, bd=bd, al=al)
    # p = drawAngleArc(p=p, r=50, aMin=0, aMax=al, aLabel=deg2rad(0), label="α", linecolor=:gray)
    plotInvolute(axs=ax, bd=bd, al=al, thMax=althMax, linecolor=:blue )
    # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(31), linecolor=:blue )
    # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(106), linecolor=:blue )
    plotInvoluteConstruction(axs=ax, bd=bd, al=al, th=althMax, linecolor=:cyan, linestyle=:solid )
    # p = drawAngleArc(p=p, r=60, aMin=0, aMax=althMax, aLabel=deg2rad(0), label=L"\theta_\alpha", linecolor=:gray)
    # savefig("on220826_involutes_alpha.svg")

    plotAlpha(axs=ax, bd=bd, al=bt)
    # p = drawAngleArc(p=p, r=30, aMin=0, aMax=bt, aLabel=deg2rad(0), label="β", linecolor=:gray)
    plotInvolute(axs=ax, bd=bd, al=bt, thMax=btthMax, linecolor=:red )
    plotInvoluteConstruction(axs=ax, bd=bd, al=bt, th=btthMax, linecolor=:magenta, linestyle=:solid )
    # p = drawAngleArc(p=p, r=70, aMin=0, aMax=btthMax, aLabel=deg2rad(0), label=L"\theta_\beta", linecolor=:gray)
    # savefig("on220826_involutes_beta.svg")

    drawToothTop(axs=ax, bd=bd, al=al, thal=althMax, bt=bt, thbt=btthMax, linecolor=:orange, linestyle=:solid)
  end



end #InvoluteTooth
