function psi_final = findFinalHeading(nextGate, afterNextGate)

vREL = afterNextGate - nextGate;
if (vREL(1) < 0 && vREL(2) > 0)||(vREL(1) < 0 && vREL(2) < 0)
    psi_final =(pi + atan(vREL(2)/vREL(1)))/2;
else
    psi_final = atan(vREL(2)/vREL(1))/2;
end

end