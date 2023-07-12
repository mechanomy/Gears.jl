module Gears
  # References
  # 2004_oberg p1599
  # 1989_Hindhede_SpurGearFundamentals - has good diagrams
  # 2012_Dooner_KinematicGeometryOfGearing

  using UnitTypes

  # structs to do conversions

  UnitTypes.@makeMeasure Pitch UnitTypes.AbstractMeasure 1.0 # units of 1/inch
  # PerInch?

  # struct GearModule
  # end

  """
  tooth involutes are perpendicular to this circle
  in meshing gears, connecting two base circles with a tangent line will cross the involute teeth at the point of contact
  """
  # baseDiameter(pitchDiameter::T where T<:AbstractDiameter, pressureAngle::U where U<: AbstractAngle) = Diameter(convert(Diameter{Inch}, pitchDiameter) * cos(pressureAngle)) :: V where V<:AbstractDiameter
  # baseDiameter(pitchDiameter::T where T<:AbstractDiameter, pressureAngle::U where U<: AbstractAngle) = Diameter(convert(Diameter{Inch}, pitchDiameter) * cos(pressureAngle)) :: V where V<:AbstractDiameter
  baseDiameter(pitchDiameter::Diameter{Inch}, pressureAngle::U where U<: AbstractAngle) = ( pitchDiameter * cos(pressureAngle))

  """
  the number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch
  """
  # diametralPitch(circularPitch) = π/circularPitch
  diametralPitch(nTeeth::Int, pitchDiameter::Diameter) = Pitch(nTeeth/convert(Float64, convert(Inch,pitchDiameter)))
  # diametralPitch(nTeethPinion, gearRatio, centerDistance) = nTeethPinion*(gearRatio+1)/2/centerDistance

  """
  the diameter of the pitch circle
  """
  pitchDiameter(dp::Pitch, nTeeth::Int) = Diameter(Inch(nTeeth / dp)) :: T where T<:AbstractDiameter # note this is assuming dp as 1/inch, would rather have this convert...
  pitchDiameter(nTeeth::Int, dp::Pitch) = Diameter(Inch(nTeeth / dp)) :: T where T<:AbstractDiameter # note this is assuming dp as 1/inch, would rather have this convert...

  nTeeth(dp::Pitch, pitchDiameter::T) where T<:AbstractDiameter = convert(Int64,dp*convert(Inch,pitchDiameter))

  # nTeeth(pitchD::T, circularPitch) = π*pitchDiameter/circularPitch

  """
  height of tooth above pitch circle
  """
  addendum(diametral::Pitch) = Inch(1/convert(Float64, diametral.value))

  """
  depth of tooth below pitch circle
  """
  dedendum(diametral::Pitch) = Inch(1.250/convert(Float64, diametral.value))


  """
  """
  outsideDiameter(pitchDiameter::UnitTypes.Diameter, addendum) = pitchDiameter + Diameter(2*addendum)


  """

  """
  rootDiameter(pitchDiameter::UnitTypes.Diameter, dedendum) = pitchDiameter - Diameter(2*dedendum)

  """
  distance between successive involutes measured along the segment of the base circle, mating teeth must have the same base pitch length
  """
  # basePitch(nTeeth::Int, baseDiameter::Diameter) = π*baseDiameter/nTeeth # length/tooth
  baseInterval(nTeeth::Int, baseDiameter::Diameter) = π*convert(typeof(baseDiameter).parameters[1],baseDiameter)/nTeeth # length/tooth
  
  # """
  #   the arc distance between similar points on the adjacent teeth in linear units along the arc segment
  # """
  # circularPitch(pitchDiameter, nTeeth) = π*pitchDiameter / nTeeth
  # # circularPitch(diametralPitch) = π/diametralPitch

  # """
  # module is the amount of pitch diameter per tooth, higher module larger tooth
  # """
  # gearModule(pitchDiameterMM, nTeeth) = pitchDiameterMM/nTeeth

  # gearRatio(nTeethPinion, nTeethGear) = nTeethGear / nTeethPinion

  abstract type AbstractGear end
  """
  Struct to represent ANSI-specification spur gears
  """
  struct GearANSI <: AbstractGear
    nTeeth::Int
    pitch::T where T<:UnitTypes.AbstractDiameter
    pressure::U where U<:UnitTypes.AbstractAngle
    diametral::Pitch
    addendum::V where V<:UnitTypes.AbstractExtent
    dedendum::W where W<:UnitTypes.AbstractExtent
    outside::X where X<:UnitTypes.AbstractDiameter
    base::Y where Y<:UnitTypes.AbstractDiameter
    root::Z where Z<:UnitTypes.AbstractDiameter
    baseInterval::A where A<:UnitTypes.AbstractExtent

    # bore::UnitTypes.AbstractDiameter
    # circularPitch::Pitch
    # wholeDepth::UnitTypes.AbstractExtent
    # workingDepth::UnitTypes.AbstractExtent
    # clearance::UnitTypes.AbstractExtent
  end

  """
  """
  function GearANSI(dp::Pitch, n::Int, pressure::UnitTypes.Degree=UnitTypes.Degree(20))#, bore::Diameter{Inch}=0 ) 
    pitch = pitchDiameter(dp, n)
    bd = baseDiameter(pitch,pressure)
    ad = addendum(dp)
    dd = dedendum(dp)
    od = outsideDiameter(pitch, ad)
    rd = rootDiameter(pitch, dd)
    bi = baseInterval(n, bd)
    return GearANSI( n, pitch, pressure, dp, ad, dd, od, bd, rd, bi)
  end
  GearANSI(n::Int, dp::Pitch, pressure::UnitTypes.Degree=UnitTypes.Degree(20)) = GearANSI(dp, n, pressure)

  function GearANSI(n::Int, pitch::T, pressure::UnitTypes.Degree=UnitTypes.Degree(20)) where T<:AbstractDiameter#, bore::Diameter{Inch}=0 ) 
    bd = baseDiameter(pitch, pressure)
    dp = diametralPitch(n, pitch)
    ad = addendum(dp)
    dd = dedendum(dp)
    od = outsideDiameter(pitch, ad)
    rd = rootDiameter(pitch, dd)
    bi = baseInterval(n, bd)
    return GearANSI( n, pitch, pressure, dp, ad, dd, od, bd, rd, bi)
  end

  function GearANSI(dp::Pitch, pitch::T, pressure::UnitTypes.Degree=UnitTypes.Degree(20)) where T<:AbstractDiameter#, bore::Diameter{Inch}=0 ) 
    # bd = baseDiameter(pitch, pressure)
    # dp = diametralPitch(n, pitch)
    # n = nTeeth(dp, pitch)
    n = round(dp.value * pitch.value)

    ad = addendum(dp)
    dd = dedendum(dp)
    od = outsideDiameter(pitch, ad)
    rd = rootDiameter(pitch, dd)
    bi = baseInterval(n, bd)
    return GearANSI( n, pitch, pressure, dp, ad, dd, od, bd, rd, bi)
  end



  """
      pulley2String(p::PlainPulley) :: String
    Returns a descriptive string of the given PlainPulley `p` of the form:
      PlainPulley[struct] @ [1.000mm,2.000mm] r[3.000mm] arrive[57.296°] depart[114.592°] aWrap[57.296°] lWrap[3.000mm]"
  """
  function gear2string(g::GearANSI)::String 
    # return @sprintf("GearANSI:\n nTeeth: %d \n pressureAngleDegree: %3.3f[°] \n pitchDiameterInch: %3.3f[in], diametralPitch: %3.3f[1/in]",
    #   g.nTeeth,
    #   convert(Float64, g.pressureAngle.value),
    #   convert(Float64, g.pitchDiameter.value),
    #   convert(Float64, g.diametralPitch.value)
    # )
    return "GearANSI: $(g.nTeeth) tooth $(g.pitch) $(g.diametral)"
  end

  """
      Base.show(io::IO, p::AbstractPulley)
    Function to `show()` a AbstractPulley via [`pulley2String`](#BeltTransmission.pulley2String).
  """
  function Base.show(io::IO, g::GearANSI)
    print(io, gear2string(g))
  end






  # nTeeth diameter to module .. 

  # plotGear() using InvoluteTooth
    # function plotTooth(; r, gm, dl, ht)

  module InvoluteTooth
    using Plots
    using Roots
    using LaTeXStrings
    using UnitTypes
    using ..Gears

    rad2deg(rad) = rad*180/π

    """
    calculates the x coordinate of a point along the involute
    Note that the handedness of the involute is determined by th-al: if th-al>0 a right hand involute, th-al<0 is a left hand
    """
    ix(r,al,th) = r*( cos(th) + (th-al)*sin(th) )

    """
    calculates the y coordinate of a point along the involute
    """
    iy(r,al,th) = r*( sin(th) - (th-al)*cos(th) )

    """
    finds thMax to intersect with the line of symmetry gamma 
    Gamma is the angle of the tooth's line of symmetry
    Alpha is the angle of the involutes root on the base circle
    Theta0 is the guess of the angle corresponding to the intersection of the involute and gamma
    """
    function calcThetaIntersect(; gm, al, th0=gm+(gm-al)*2)
      #calculate thMax to intersect with the gamma line
      fxi(thc,al,gm) = cos(thc)^2 + 2*(thc-al)*cos(thc)*sin(thc) + (sin(thc)^2-cos(gm)^2)*(thc-al)^2 - cos(gm)^2
      fyi(thc,al,gm) = sin(thc)^2 - 2*(thc-al)*cos(thc)*sin(thc) + (cos(thc)^2-sin(gm)^2)*(thc-al)^2 - sin(gm)^2
      fxz(thc) = fxi(thc, al, gm)^2 + fyi(thc, al, gm)^2 
      return find_zero(fxz, th0, atol=1e-5) 
    end

   function calcThetaRootIntersection(; rrd, rbd, alb, al0 )
      # fbt(thc) = rrd - rbd*sqrt(1+(thc-btb)^2)
      # thbtr = find_zero(fbt, bto, atol=1e-5) 

      fal(thc) = rrd - rbd*sqrt(1+(thc-alb)^2)
      # return find_zero(fal, alb, atol=1e-5) 
      return find_zero(fal,al0, atol=1e-5) 
    end

    """
    Finds the angle when the involute has risen through toothHeight
    Gamma is the angle of the tooth's line of symmetry
    Alpha is the angle of the involute's root on the base circle
    toothHeight is the desired radial height of the tooth beyond the base circle
    """
    function calcThetaHeight(; r, gm, al, toothHeight)
      htal(th) = r * ( cos(th-gm) + (th-al)*sin(th-gm) - 1 )
      fal(th) = toothHeight - htal(th)
      thal = find_zero(fal, gm + (gm-al)*2)
      # @show rad2deg(thal)
      return thal
    end

    """
    Plots the base arc under the tooth
    """
    function plotBaseSection(; r, gm, al, bt, thExtra=0.5, n=100, p=nothing, linecolor=:gray, linestyle=:dash)
      if p == nothing
        p = plot(aspect_ratio=:equal)
      else
        p = plot(p, aspect_ratio=:equal)
      end

      #want to plot only a section of the base circle, that from al to thMax
      ths = LinRange(al-thExtra, bt+thExtra, n) 
      xs = r.*cos.(ths)
      ys = r.*sin.(ths)
      p = plot(p, xs, ys, linecolor=linecolor, linestyle=linestyle)
      return p
    end

    """
    Plots the line of symmetry gamma
    """
    function plotGamma(; p, gm, r, ht, linecolor=:green, linestyle=:dash )
      p = plot(p, [0,(r+ht)*cos(gm)], [0,(r+ht)*sin(gm)], linecolor=linecolor, linestyle=linestyle)
      return p
    end

    """
    Plots the line to the base of the involute
    """
    function plotAlpha(; p, r, al , linecolor=:gray, linestyle=:dash)
      p = plot(p, [0,r*cos(al)], [0,r*sin(al)], linecolor=linecolor, linestyle=linestyle)
      return p
    end

    """
    Plots an involute of the circle
    """
    function plotInvolute(; r, gm, al, thMax, p=nothing, linecolor=:blue, linestyle=:solid)
      if p === nothing
        p = plot(aspect_ratio=:equal)
      else
        p = plot(p, aspect_ratio=:equal)
      end

      ths = LinRange(al, thMax,100) 

      xs = ix.(r, al, ths)
      ys = iy.(r, al, ths)
      p = plot!(xs,ys, linecolor=linecolor, linestyle=linestyle)

      p = plot(p, legend=false)
      return p
    end

    """
    """
    function plotInvoluteConstruction(; r, al, th, p=nothing, linecolor=:gray, linestyle=:dash)
      xs = [0, r*cos(th), ix(r,al,th)]
      ys = [0, r*sin(th), iy(r,al,th)]
      p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle)
    end

    """
    """
    function drawAngleArc(; p, r, aMax, aMin=0, label="label", aLabel=aMax/2, linecolor=:gray, linestyle=:solid, fontsize=10)
      ths = LinRange(aMin, aMax, 100)
      xs = r .* cos.(ths)
      ys = r .* sin.(ths)
      p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle, markerstyle=:rightarrow)

      # annotate!(r*cos(aLabel), r*sin(aLabel), "  "*label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationrotation=rad2deg(aLabel) )
      annotate!(r*cos(aLabel), r*sin(aLabel), label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationhalign=:left, annotationvalign=:bottom )
      return p
    end

    """
    """
    function drawToothTop(;r, al, thal, bt, thbt, p=nothing, linecolor=:orange, linestyle=:dash)
      p = plot(p, [ix(r,al,thal), ix(r,bt,thbt)], [iy(r,al,thal), iy(r,bt,thbt)], linecolor=linecolor, linestyle=linestyle)
      return p
    end

    """
    plots both sides to form a tooth
    r is the root diameter
    gamma is the center angle of the tooth
    delta is the angular width of the tooth at the root diameter
    height is the tooth height, or the difference between the outside and root diameters
    """
    function plotTooth(; r, gm, dl, ht, p=nothing)
      # plot construction of both involutes
      al = gm - dl/2
      althMax = calcThetaHeight(r=r, gm=gm, al=al, ht=ht)
      bt = gm + dl/2
      btthMax = calcThetaHeight(r=r, gm=gm, al=bt, ht=ht)

      if p === nothing
        p = plot(legend=false, title="Involute tooth with γ=$(round(rad2deg(gm))), δ=$(round(rad2deg(dl)))")
      end
      p = plot(p, [0,r*1.1],[0,0], linecolor=:black)
      p = plotBaseSection(p=p, r=r, gm=gm, al=al, bt=bt, linecolor=:black, linestyle=:solid)
      p = plotGamma(p=p, r=r, gm=gm, ht=ht)
      # p = drawAngleArc(p=p, r=40, aMin=0, aMax=gm, aLabel=deg2rad(0), label="γ", linecolor=:gray)
      # savefig("on220826_involutes_baseGamma.svg")

      # p = drawAngleArc(p=p, r=70, aMin=al, aMax=bt, aLabel=deg2rad(35), label=L"\delta", linecolor=:gray)
      # savefig("on220826_involutes_delta.svg")

      p = plotAlpha(p=p, r=r, al=al)
      # p = drawAngleArc(p=p, r=50, aMin=0, aMax=al, aLabel=deg2rad(0), label="α", linecolor=:gray)
      p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=althMax, linecolor=:blue )
      # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(31), linecolor=:blue )
      # p = plotInvolute(p=p, r=r, gm=gm, al=al, thMax=deg2rad(106), linecolor=:blue )
      p = plotInvoluteConstruction(p=p, r=r, al=al, th=althMax, linecolor=:cyan, linestyle=:solid )
      # p = drawAngleArc(p=p, r=60, aMin=0, aMax=althMax, aLabel=deg2rad(0), label=L"\theta_\alpha", linecolor=:gray)
      # savefig("on220826_involutes_alpha.svg")

      p = plotAlpha(p=p, r=r, al=bt)
      # p = drawAngleArc(p=p, r=30, aMin=0, aMax=bt, aLabel=deg2rad(0), label="β", linecolor=:gray)
      p = plotInvolute(p=p, r=r, gm=gm, al=bt, thMax=btthMax, linecolor=:red )
      p = plotInvoluteConstruction(p=p, r=r, al=bt, th=btthMax, linecolor=:magenta, linestyle=:solid )
      # p = drawAngleArc(p=p, r=70, aMin=0, aMax=btthMax, aLabel=deg2rad(0), label=L"\theta_\beta", linecolor=:gray)
      # savefig("on220826_involutes_beta.svg")

      p = drawToothTop(p=p, r=r, al=al, thal=althMax, bt=bt, thbt=btthMax, linecolor=:orange, linestyle=:solid)

    return p
    end

    function getToothProfilePoints(g::G where G<:Gears.AbstractGear; nPerTooth::Int=100)
      nPerInvolute = convert(Int64,round(nPerTooth/2)) #this may be better defined by physical spacing...

      rod = convert(Float64, convert(Inch, g.outside/2))
      rpd = convert(Float64, convert(Inch, g.pitch/2))
      rrd = convert(Float64, convert(Inch, g.root/2))
      rbd = convert(Float64, convert(Inch, g.base/2))
      t = π/2/24
    
      xs = []
      ys = []
      for i in 1:convert(Int64, round(g.nTeeth))
        gm=2*π/g.nTeeth*(i) #gamma is the angle of the line of the tooth center

        psi = acos(rbd/rpd)
        vpsi = tan(psi)-psi

        alb = gm - asin(t/2/rpd) - vpsi
        alo = Gears.InvoluteTooth.calcThetaHeight(r=rbd, gm=gm, al=alb, toothHeight=rod-rbd) # this is the angle of the involute tip at the outside diameter
        thalr = Gears.InvoluteTooth.calcThetaRootIntersection(rrd=rrd, rbd=rbd, alb=alb, al0=alb)
        ths = LinRange(thalr, alo, nPerInvolute)
        append!(xs, Gears.InvoluteTooth.ix.(rbd, alb, ths) )
        append!(ys, Gears.InvoluteTooth.iy.(rbd, alb, ths) )

        btb = gm + asin(t/2/rpd) + vpsi
        bto = Gears.InvoluteTooth.calcThetaHeight(r=rbd, gm=gm, al=btb, toothHeight=rod-rbd) 
        thbtr = Gears.InvoluteTooth.calcThetaRootIntersection(rrd=rrd, rbd=rbd, alb=btb, al0=bto)
        ths = LinRange(bto, thbtr, nPerInvolute)  # reverse direction
        append!(xs, Gears.InvoluteTooth.ix.(rbd, btb, ths) )
        append!(ys, Gears.InvoluteTooth.iy.(rbd, btb, ths) )

        # alr = convert(Float64,convert(Radian,gm)) - (π/2/g.nTeeth + tan(convert(Float64,convert(Radian,g.pressure)))-convert(Float64,convert(Radian,g.pressure))) #alpha is the angle of the root of the involute preceeding gamma
        # alo = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=alr, toothHeight=rod-rrd) # this is the angle of the involute tip at the outside diameter
        # # p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=al, thMax=althPdOd, p=p, linecolor=:blue, linestyle=:solid)
        # ths = LinRange(alr, alo, nPerInvolute) 
        # append!(xs, Gears.InvoluteTooth.ix.(rrd, alr, ths) )
        # append!(ys, Gears.InvoluteTooth.iy.(rrd, alr, ths) )

        # btr = convert(Float64,convert(Radian,gm)) + (π/2/g.nTeeth + tan(convert(Float64,convert(Radian,g.pressure)))-convert(Float64,convert(Radian,g.pressure))) #beta is the angle of the involute suceeding gamma
        # bto = Gears.InvoluteTooth.calcThetaHeight(r=rrd, gm=gm, al=btr, toothHeight=rod-rrd ) 
        # # p = Gears.InvoluteTooth.plotInvolute( r=rrd, gm=gm, al=bt, thMax=btthPdOd, p=p, linecolor=:red, linestyle=:solid)
        # ths = LinRange(bto, btr, nPerInvolute)  # reverse direction
        # append!(xs, Gears.InvoluteTooth.ix.(rrd, btr, ths) )
        # append!(ys, Gears.InvoluteTooth.iy.(rrd, btr, ths) )
      end
      append!(xs, xs[1])
      append!(ys, ys[1])

      # p = plot(xs, ys, linecolor=:blue, linestyle=:solid, markerstyle=:square, markersize=3, markercolor=:red, aspect_ratio=:equal)
      # display(p)
      
      return (xs,ys)
    end

    function writeToothProfilePoints(g::G where G<:Gears.AbstractGear; fileName::String="", fileExtension::String="txt", nPerTooth::Int=100)
      (xs, ys) = getToothProfilePoints(g, nPerTooth=nPerTooth)

      #now emit the points in the format expected by CAD
      zs = xs *0

      if fileName == ""
        fileName = "gearProfilePoints_$(g.diametral)_$(g.nTeeth).sldcrv"
      end

      #solidworks format: 25mmTAB0mmTAB0mmCRLF
      open(fileName, "w") do f
        for i in 1:length(xs)
          # write(f, @sprintf("%3.3fmm\t%3.3fmm\t%3.3fmm", xs[i],ys[i],zs[i]))
          # print("$(xs[i])mm\t$(ys[i])mm\t$(zs[i])mm\r\n")
          # write(f, "$(xs[i])mm\t$(ys[i])mm\t$(zs[i])mm\r\n")
          write(f, "$(xs[i])in\t$(ys[i])in\t$(zs[i])in\r\n")
        end
      end
    end

  end #InvoluteTooth
end #Gears

