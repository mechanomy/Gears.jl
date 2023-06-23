# module Gears
#   struct ExternalSpurGear #2004_oberg p1599
#     diametralPitch::Length
#   end
# end

module InvoluteTooth
using Plots

rad2deg(rad) = rad*180/π

ix(r,al,th) = r*( cos(th) + (th-al)*sin(th) )
iy(r,al,th) = r*( sin(th) - (th-al)*cos(th) )

function calcThetaIntersect(;gm, al)
  return gm + (gm-al)*2
end

using Roots
function calcThetaIntersect(; gm, al, th0=gm+(gm-al)*2)
  #calculate thMax to intersect with the gamma line
  fxi(thc,al,gm) = cos(thc)^2 + 2*(thc-al)*cos(thc)*sin(thc) + (sin(thc)^2-cos(gm)^2)*(thc-al)^2 - cos(gm)^2
  fyi(thc,al,gm) = sin(thc)^2 - 2*(thc-al)*cos(thc)*sin(thc) + (cos(thc)^2-sin(gm)^2)*(thc-al)^2 - sin(gm)^2
  fxz(thc) = fxi(thc, al, gm)^2 + fyi(thc, al, gm)^2 
  return find_zero(fxz, th0, atol=1e-5) 
end

function calcThetaHeight(; r, gm, al, ht)
  htal(th) = r * ( cos(th-gm) + (th-al)*sin(th-gm) - 1 )
  fal(th) = ht - htal(th)
  # thal = find_zero(fal, (gm+al)/2)
  thal = find_zero(fal, gm + (gm-al)*2)
  # @show rad2deg(thal)
  return thal
end


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

function plotGamma(; p, gm, r, ht, linecolor=:green, linestyle=:dash )
  p = plot(p, [0,(r+ht)*cos(gm)], [0,(r+ht)*sin(gm)], linecolor=linecolor, linestyle=linestyle)
  return p
end
function plotAlpha(; p, r, al , linecolor=:gray, linestyle=:dash)
  p = plot(p, [0,r*cos(al)], [0,r*sin(al)], linecolor=linecolor, linestyle=linestyle)
  return p
end

"""
Plots an involute of the circle
"""
function plotInvolute(; r, gm, al, thMax, p=nothing, linecolor=:blue, linestyle=:solid)
  if p == nothing
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

function plotInvoluteConstruction(; r, al, th, p=nothing, linecolor=:gray, linestyle=:dash)
  xs = [0, r*cos(th), ix(r,al,th)]
  ys = [0, r*sin(th), iy(r,al,th)]
  p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle)
end

using LaTeXStrings
function drawAngleArc(; p, r, aMax, aMin=0, label="label", aLabel=aMax/2, linecolor=:gray, linestyle=:solid, fontsize=10)
  ths = LinRange(aMin, aMax, 100)
  xs = r .* cos.(ths)
  ys = r .* sin.(ths)
  p = plot(p, xs,ys, linecolor=linecolor, linestyle=linestyle, markerstyle=:rightarrow)

  # annotate!(r*cos(aLabel), r*sin(aLabel), "  "*label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationrotation=rad2deg(aLabel) )
  annotate!(r*cos(aLabel), r*sin(aLabel), label, annotationcolor=linecolor, annotationfontsize=fontsize, annotationhalign=:left, annotationvalign=:bottom )
  return p
end

function drawToothTop(;r, al, thal, bt, thbt, p=nothing, linecolor=:orange, linestyle=:dash)
  p = plot(p, [ix(r,al,thal), ix(r,bt,thbt)], [iy(r,al,thal), iy(r,bt,thbt)], linecolor=linecolor, linestyle=linestyle)
  return p
end

"""plots both sides to form a tooth"""
function plotTooth(; r, gm, dl, ht)
  # plot construction of both involutes
  al = gm - dl/2
  # althMax = calcThetaIntersect(gm=gm, al=al)
  althMax = calcThetaHeight(r=r, gm=gm, al=al, ht=ht)
  bt = gm + dl/2
  # btthMax = calcThetaIntersect(gm=gm, al=bt)
  btthMax = calcThetaHeight(r=r, gm=gm, al=bt, ht=ht)

  p = plot(legend=false, title="Involute tooth with γ=$(round(rad2deg(gm))), δ=$(round(rad2deg(dl)))")
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

using Plots
# p = InvoluteTooth.plotTooth( r=96, gm=deg2rad(45), dl=deg2rad(56), ht=20 )
# savefig("on220826_involutes_toothHeight.svg")
# display(p)
# p = plotTooth( r=96, gm=deg2rad(45), dl=deg2rad(3) )
# display(p)
 


# # given: khk_KHK-US_PSA2-100J30
# gearmodule=2
# pressureAngle=20 #deg
# nTeeth = 100
# bore = 15
# pitchd = 200
# outsided = 204
# facewidth=20
# allowabletorque=34.2 #Nm- bending strength
# # backlash=0-0.46
# weight=0.72


# # module = pitchd / nTeeth
# addendum = gearmodule
# dedendum = 1.25*gearmodule
# pitchd = gearmodule * nTeeth
# # outsided = pitchd + 2*gearmodule = gearmodule*(nTeeth+2)
# rootd = pitchd - 2.5*gearmodule
# based = pitchd * cosd(pressureAngle)
# basepitch = gearmodule * π * cosd(pressureAngle)
# toothThick = gearmodule * π /2 #tooth thickness at standard pitch diameter
# centerdistance = m*(n1+n2)/2

#using rootd
# dl = atan(toothThick/2, rootd/2) * 2
# @show InvoluteTooth.rad2deg(dl)
# @show toothHeight = (outsided-rootd)/2
# p = InvoluteTooth.plotTooth( r=rootd/2, gm=deg2rad(90), dl=dl, ht=toothHeight )
# display(p)

# # using pitchd
# dl = atan(toothThick/2, pitchd/2) * 2
# @show InvoluteTooth.rad2deg(dl)
# @show toothHeight = (outsided-pitchd)/2
# p = InvoluteTooth.plotTooth( r=pitchd/2, gm=deg2rad(90), dl=dl, ht=toothHeight )
# display(p)

# using based
# dl = atan(toothThick/2, based/2) * 2
# InvoluteTooth.rad2deg(dl)
# toothHeight = (outsided-based)/2
# p = InvoluteTooth.plotTooth( r=based/2, gm=deg2rad(90), dl=dl, ht=toothHeight )
# display(p)

# function printEquationCurve( gearmodule, nTeeth, gm=π/2, pressureAngle=deg2rad(20) )
#   pitch = 1/gearmodule #1/mm
#   addendum = gearmodule
#   dedendum = 1.25*gearmodule
#   @show pitchD = gearmodule * nTeeth
#   @show outsideD = pitchD + 2*gearmodule # = gearmodule*(nTeeth+2)
#   @show rootD = pitchD - 2.5*gearmodule
#   @show baseD = pitchD * cos(pressureAngle)
#   @show circularpitch = gearmodule * π
#   @show basepitch = gearmodule * π * cos(pressureAngle) #pitch between gear teeth along the base circle
#   # @show toothThick = gearmodule * π /2 #tooth thickness at standard pitch diameter. this from sdpsi is /2 what I see in khk module2 gears. is defined on 2000_sdpsi_gears #15

#   # r = (pitchD - 2*dedendum)/2
#   r = baseD/2
#   # dl = atan(toothThick/2, r) * 2
#   # dl = atan(toothThick/2, pitchD/2) * 2

#   @show toothThick = π / 2 / pitch
#   # dl = 2*asin( toothThick / 2/r)
#   dl = deg2rad(3.7*2)
#   al = gm - dl/2

#   # dl = 2*π/nTeeth # this ignores the groove
#   @show rad2deg(dl)
#   # toothHeight = (outsideD-baseD)/2
#   @show toothHeight = addendum + dedendum
#   # p = InvoluteTooth.plotTooth( r=r, gm=gm, dl=dl, ht=toothHeight )
#   # display(p)

#   thal = InvoluteTooth.calcThetaHeight(r=r, gm=gm, al=al, ht=toothHeight)

#   xt = "$r*cos(t) + $r*(t-$al)*sin(t)"
#   yt = "$r*sin(t) - $r*(t-$al)*cos(t)"
#   t1 = "$al"
#   t2 = "$thal"

#   println("Equation Curve = ")
#   println(xt)
#   println(yt)
#   println(t1)
#   println(t2)
# end

function printEquationCurve( gearmodule, nTeeth, gamma=π/2, pressureAngle=deg2rad(20) )
  # pitch = 1/gearmodule #1/mm
  addendum = gearmodule
  dedendum = 1.25*gearmodule
  pitchD = gearmodule * nTeeth
  # outsideD = pitchD + 2*gearmodule 
  # rootD = pitchD - 2.5*gearmodule
  # baseD = pitchD * cos(pressureAngle)
  # circularpitch = gearmodule * π
  # basepitch = gearmodule * π * cos(pressureAngle) #pitch between gear teeth along the base circle

  # r = baseD/2 
  r = (pitchD - 2*dedendum)/2
  toothThick = gearmodule * pi/2
  # delta = deg2rad(3.7*2) 
  delta = 2*asin( toothThick / 2/r)
  alpha = gamma - delta/2

  toothHeight = addendum + dedendum
  p = InvoluteTooth.plotTooth( r=r, gm=gamma, dl=delta, ht=toothHeight )
  display(p)

  thal = InvoluteTooth.calcThetaHeight(r=r, gm=gamma, al=alpha, ht=toothHeight)

  println("Equation Curve with r=$r and δ=$delta for a $nTeeth tooth module $gearmodule gear: ")
  println("xt = $r*cos(t) + $r*(t-$alpha)*sin(t)")
  println("yt = $r*sin(t) - $r*(t-$alpha)*cos(t)")
  println("t1 = $alpha")
  println("t2 = $thal")
end



# printEquationCurve(nTeeth, based, outsided )

## sdpsi_khk_KHK-US_NSU2-32J20
nTeeth=32
gearmodule=2
pressureAngle=deg2rad(20) #deg
printEquationCurve( gearmodule, nTeeth, π/2, pressureAngle)
# @show dl = 2*π/nTeeth
# @show toothHeight = (outsided-based)/2
# p = InvoluteTooth.plotTooth( r=based/2, gm=deg2rad(90), dl=dl, ht=toothHeight )
# display(p)