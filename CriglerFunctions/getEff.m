function [eta] = getEff(Pc_k, ek, interpMode)

crieglerPlots
if strcmp(interpMode, 'linear')
    etaVec = interp1([1 0], plots.eta.etaMat', ek, 'linear', 'extrap');
    eta = interp1(plots.eta.PckVec, etaVec, Pc_k, 'linear', 'extrap');
end

if strcmp(interpMode, 'spline')
    etaVec = interp1([1 0], plots.eta.etaMat', ek, 'linear', 'extrap');
    eta = interp1(plots.eta.PckVec, etaVec, Pc_k, 'spline');
end