function [h, xref] = equalCon(U, info)
xf = info.xf(1); yf = info.xf(2); zf = info.xf(3);
psif = info.psif;
x0 = info.x0;
info.Ts = U(end);

xref = predictStates(x0,U(1:info.Nsteps*info.dimM),info);

[xi, yi, zi] = extractPos(xref, info.Nsteps+1);
psi = extractHeading(xref, info.Nsteps+1);

if ~isempty(psif)
    h = [xi(end) - xf; yi(end)-yf; zi(end)-zf; psi(end)-psif];
else
     h = [xi(end) - xf; yi(end)-yf; zi(end)-zf];
end