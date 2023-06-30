module Gears

  # References
  # 2004_oberg p1599
  # 1989_Hindhede_SpurGearFundamentals - has good diagrams
  # 2012_Dooner_KinematicGeometryOfGearing

  using UnitTypes

  # structs to do conversions

  UnitTypes.@makeMeasure Pitch UnitTypes.AbstractMeasure 1.0 # units of 1/inch

  struct GearModule
  end

  """
  tooth involutes are perpendicular to this circle
  in meshing gears, connecting two base circles with a tangent line will cross the involute teeth at the point of contact
  """
  baseCircle(pitchDiameter::UnitTypes.Diameter, pressureAngle::T where T<: UnitTypes.AbstractAngle) = pitchDiameter * cos(pressureAngle)
  """
  meshed gears have these circles tangent
  """
  # pitchCircle

  """
  distance between successive involutes measured along the segment of the base circle, mating teeth must have the same base pitch length
  """
  basePitch(nTeeth, baseDiameter) = π*baseDiameter/nTeeth
  
  """
    the arc distance between similar points on the adjacent teeth in linear units along the arc segment
  """
  circularPitch(pitchDiameter, nTeeth) = π*pitchDiameter / nTeeth
  circularPitch(diametralPitch) = π/diametralPitch

  """
  module is the amount of pitch diameter per tooth, higher module larger tooth
  """
  gearModule(pitchDiameterMM, nTeeth) = pitchDiameterMM/nTeeth

  """
  the number of teeth per inch of pitch diameter, larger DP is smaller teeth; dp is not a pitch in the same sense as basePitch
  """
  diametralPitch(circularPitch) = π/circularPitch
  diametralPitch(nTeeth,pitchDiameter) = nTeeth/pitchDiameter
  diametralPitch(nTeethPinion, gearRatio, centerDistance) = nTeethPinion*(gearRatio+1)/2/centerDistance


  centerDistance(nTeethPinion, gearRatio, diametralPitch) = nTeethPinion * (gearRatio+1) / 2 / diametralPitch
  centerDistance(pitchDiameterPinion, pitchDiameterGear) = (pitchDiameterPinion+pitchDiameterGear)/2
  centerDistance(nTeethPinion, nTeethGear, diametralPitch) = (nTeethPinion+nTeethGear)/ 2/ diametralPitch
  centerDistance(nTeethPinion, nTeethGear, circularPitch) = (nTeethPinion+nTeethGear)*circularPitch /2 /π



  gearRatio(nTeethPinion, nTeethGear) = nTeethGear / nTeethPinion

  nTeeth(diametralPitch, pitchDiameter) = diametralPitch*pitchDiameter
  nTeeth(pitchDiameter,circularPitch) = π*pitchDiameter/circularPitch

  """
  height of tooth above pitch circle
  """
  addendum(diametral::Pitch) = 1/diametral
  """
  depth of tooth below pitch circle
  """
  function dedendum(diametral::Pitch) :: T where T<:UnitTypes.AbstractExtent
    return 1/diametral
  end

  """
  Struct to represent ANSI-specification spur gears
  """
  struct SpurANSI
    nTeeth::Int
    pressureAngle::UnitTypes.AbstractAngle
    diametralPitch::Pitch
    # pitchDiameter::UnitTypes.Inch
  end
  # struct SpurANSI
  #   nTeeth::Int
  #   pitch::UnitTypes.AbstractDiameter
  #   base::UnitTypes.AbstractDiameter
  #   root::UnitTypes.AbstractDiameter
  #   outside::UnitTypes.AbstractDiameter
  #   bore::UnitTypes.AbstractDiameter
  #   circularPitch::Pitch
  #   diametralPitch::UnitTypes.AbstractExtent
  #   pressureAngle::UnitTypes.AbstractAngle
  #   wholeDepth::UnitTypes.AbstractExtent
  #   workingDepth::UnitTypes.AbstractExtent
  #   addendum::UnitTypes.AbstractExtent
  #   dedendum::UnitTypes.AbstractExtent
  #   clearance::UnitTypes.AbstractExtent
  # end

  # nTeeth diameter to module .. 

  # plotGear() using InvoluteTooth
    # function plotTooth(; r, gm, dl, ht)

  module InvoluteTooth
    using Plots
    using Roots
    using LaTeXStrings
    using UnitTypes

    rad2deg(rad) = rad*180/π

    """
    calculates the x coordinate of a point along the involute
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
    # needs Base.iterate()? https://chifi.dev/how-to-build-your-own-iterate-in-julia-56c47541879c
    # function calcThetaIntersect(; gm::T, al::T, th0::T=gm+(gm-al)*2) where {T<:UnitTypes.AbstractAngle}
    #   #calculate thMax to intersect with the gamma line
    #   fxi(thc,al,gm) = cos(thc)^2 + 2*(thc-al)*cos(thc)*sin(thc) + (sin(thc)^2-cos(gm)^2)*(thc-al)^2 - cos(gm)^2
    #   fyi(thc,al,gm) = sin(thc)^2 - 2*(thc-al)*cos(thc)*sin(thc) + (cos(thc)^2-sin(gm)^2)*(thc-al)^2 - sin(gm)^2
    #   fxz(thc) = fxi(thc, al, gm)^2 + fyi(thc, al, gm)^2 
    #   return find_zero(fxz, th0, atol=1e-5) 
    # end

    """
    finds the angle when the involute has risen to height
    Gamma is the angle of the tooth's line of symmetry
    Alpha is the angle of the involutes root on the base circle
    Height is the desired radial height of the tooth beyond the base circle
    """
    function calcThetaHeight(; r, gm, al, ht)
      htal(th) = r * ( cos(th-gm) + (th-al)*sin(th-gm) - 1 )
      fal(th) = ht - htal(th)
      # thal = find_zero(fal, (gm+al)/2)
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

  end #InvoluteTooth
end #Gears

