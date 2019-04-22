function [h, xref] = equalConFixedTime(U, info)
xf = info.xf(1); yf = info.xf(2); zf = info.xf(3);
psif = info.psif;
x0 = info.x0;

xref = predictStates(x0,U,info);

[xi, yi, zi] = extractPos(xref, info.Nsteps+1);
psi = extractHeading(xref, info.Nsteps+1);

if ~isempty(psif)
    h = [xi(end) - xf; yi(end)-yf; zi(end)-zf; psi(end)-psif];
else
     h = [xi(end) - xf; yi(end)-yf; zi(end)-zf];
end