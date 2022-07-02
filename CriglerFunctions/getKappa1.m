function [k] = getKappa1(nP, interpMode, B, Jw)

crieglerPlots

if strcmp(interpMode,'linear')
    if strcmp(nP,'singleRotation')
        if B == 2
            k = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.kVec,Jw,'linear', 'extrap');
        elseif B == 3
            k = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.kVec,Jw,'linear', 'extrap');
        elseif B == 4
            k = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.kVec,Jw,'linear', 'extrap');
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            k = interp1(plots.dualRotation.B4.JkVec,plots.dualRotation.B4.kVec,Jw,'linear', 'extrap');
        elseif B == 8
            k = interp1(plots.dualRotation.B8.JkVec,plots.dualRotation.B8.kVec,Jw,'linear', 'extrap');
        elseif B == 12
            k = interp1(plots.dualRotation.B12.JkVec,plots.dualRotation.B12.kVec,Jw,'linear', 'extrap');
        end
    end
end


if strcmp(interpMode,'spline')
    if strcmp(nP,'singleRotation')
        if B == 2
            k = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.kVec,Jw,'spline');
        elseif B == 3
            k = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.kVec,Jw,'spline');
        elseif B == 4
            k = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.kVec,Jw,'spline');
            
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            k = interp1(plots.dualRotation.B4.JkVec,plots.dualRotation.B4.kVec,Jw,'spline');
        elseif B == 8
            k = interp1(plots.dualRotation.B8.JkVec,plots.dualRotation.B8.kVec,Jw,'spline');
        elseif B == 12
            k = interp1(plots.dualRotation.B12.JkVec,plots.dualRotation.B12.kVec,Jw,'spline');
        end
    end
end

