function [k,ek] = getKappa(nP, interpMode, B, Jw)

crieglerPlots

if strcmp(interpMode,'linear')
    if strcmp(nP,'singleRotation')
        if B == 2
            k = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.ekVec,Jw,'linear', 'extrap');
        elseif B == 3
            k = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.ekVec,Jw,'linear', 'extrap');
        elseif B == 4
            k = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.ekVec,Jw,'linear', 'extrap');
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            k = interp1(plots.dualRotation.B4.JkVec,plots.dualRotation.B4.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.dualRotation.B4.JekVec,plots.dualRotation.B4.ekVec,Jw,'linear', 'extrap');
        elseif B == 8
            k = interp1(plots.dualRotation.B8.JkVec,plots.dualRotation.B8.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.dualRotation.B8.JekVec,plots.dualRotation.B.ekVec,Jw,'linear', 'extrap');
        elseif B == 12
            k = interp1(plots.dualRotation.B12.JkVec,plots.dualRotation.B12.kVec,Jw,'linear', 'extrap');
            ek = interp1(plots.dualRotation.B12.JekVec,plots.dualRotation.B12.ekVec,Jw,'linear', 'extrap');
        end
    end
end


if strcmp(interpMode,'spline')
    if strcmp(nP,'singleRotation')
        if B == 2
            k = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.kVec,Jw,'spline');
            ek = interp1(plots.singleRotation.B2.JwVec,plots.singleRotation.B2.ekVec,Jw,'spline');
        elseif B == 3
            k = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.kVec,Jw,'spline');
            ek = interp1(plots.singleRotation.B3.JwVec,plots.singleRotation.B3.ekVec,Jw,'spline');
        elseif B == 4
            k = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.kVec,Jw,'spline');
            ek = interp1(plots.singleRotation.B4.JwVec,plots.singleRotation.B4.ekVec,Jw,'spline');
            
        end
    elseif strcmp(nP,'dualRotation')
        if B == 4
            k = interp1(plots.dualRotation.B4.JkVec,plots.dualRotation.B4.kVec,Jw,'spline');
            ek = interp1(plots.dualRotation.B4.JekVec,plots.dualRotation.B4.ekVec,Jw,'spline');
        elseif B == 8
            k = interp1(plots.dualRotation.B8.JkVec,plots.dualRotation.B8.kVec,Jw,'spline');
            ek = interp1(plots.dualRotation.B8.JekVec,plots.dualRotation.B8.ekVec,Jw,'spline');
        elseif B == 12
            k = interp1(plots.dualRotation.B12.JkVec,plots.dualRotation.B12.kVec,Jw,'spline');
            ek = interp1(plots.dualRotation.B12.JekVec,plots.dualRotation.B12.ekVec,Jw,'spline');
        end
    end
end

