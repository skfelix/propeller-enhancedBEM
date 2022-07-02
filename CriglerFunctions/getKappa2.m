function [ek] = getKappa2(nP, interpMode, B, Jw)

crieglerPlots

if strcmp(interpMode,'linear')
    if strcmp(nP,'singleRotation')
        if B == 2
            ek = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.ekVec,Jw,'linear', 'extrap');
        elseif B == 3
            ek = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.ekVec,Jw,'linear', 'extrap');
        elseif B == 4
            ek = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.ekVec,Jw,'linear', 'extrap');
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            ek = interp1(plots.dualRotation.B4.JekVec,plots.dualRotation.B4.ekVec,Jw,'linear', 'extrap');
        elseif B == 8
            ek = interp1(plots.dualRotation.B8.JekVec,plots.dualRotation.B8.ekVec,Jw,'linear', 'extrap');
        elseif B == 12
            ek = interp1(plots.dualRotation.B12.JekVec,plots.dualRotation.B12.ekVec,Jw,'linear', 'extrap');
        end
    end
end


if strcmp(interpMode,'spline')
    if strcmp(nP,'singleRotation')
        if B == 2
            ek = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.ekVec,Jw,'spline');
        elseif B == 3
            ek = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.ekVec,Jw,'spline');
        elseif B == 4
            ek = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.ekVec,Jw,'spline');
            
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            ek = interp1(plots.dualRotation.B4.JekVec,plots.dualRotation.B4.ekVec,Jw,'spline');
        elseif B == 8
            ek = interp1(plots.dualRotation.B8.JekVec,plots.dualRotation.B8.ekVec,Jw,'spline');
        elseif B == 12
            ek = interp1(plots.dualRotation.B12.JekVec,plots.dualRotation.B12.ekVec,Jw,'spline');
        end
    end
end

