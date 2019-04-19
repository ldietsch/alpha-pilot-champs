function [uref, xref] = getNextTraj(Xref, Uref, info)

uref = Uref(info.currentGate+1:info.nextGateID*info.dimM*info.nMPC,1);
xref = Xref(info.currentGate+1:info.nextGateID*3*info.nMPC,1);

end
