function [KxVec] = getKx(nP, interpMode, B, Jw)

crieglerPlots

if strcmp(interpMode,'linear')
    if strcmp(nP,'singleRotation')
        if B == 2
            JwVec = plots.singleRotation.B2.JKxVec;
            KxMat = plots.singleRotation.B2.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        elseif B == 3
            JwVec = plots.singleRotation.B3.JKxVec;
            KxMat = plots.singleRotation.B3.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        elseif B == 4
            JwVec = plots.singleRotation.B4.JKxVec;
            KxMat = plots.singleRotation.B4.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        end
        
    elseif strcmp(nP,'dualRotation')
        if B == 4
            JwVec = plots.dualRotation.B4.JKxVec;
            KxMat = plots.dualRotation.B4.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        elseif B == 8
            JwVec = plots.dualRotation.B8.JKxVec;
            KxMat = plots.dualRotation.B8.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        elseif B == 12
            JwVec = plots.dualRotation.B12.JKxVec;
            KxMat = plots.dualRotation.B12.KxMat;
            KxVec = interp1(JwVec,KxMat, Jw, 'linear', 'extrap');
        end
    end
end

% spline

if strcmp(interpMode,'spline')
    if strcmp(nP,'singleRotation')
        if B == 2
            JwVec = plots.singleRotation.B2.JKxVec;
            KxMat = plots.singleRotation.B2.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        elseif B == 3
            JwVec = plots.singleRotation.B3.JKxVec;
            KxMat = plots.singleRotation.B3.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        elseif B == 4
            JwVec = plots.singleRotation.B4.JKxVec;
            KxMat = plots.singleRotation.B4.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        end
        
    elseif strcmp(nP,'dualRotation')
        if B == 4
            JwVec = plots.dualRotation.B4.JKxVec;
            KxMat = plots.dualRotation.B4.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        elseif B == 8
            JwVec = plots.dualRotation.B8.JKxVec;
            KxMat = plots.dualRotation.B8.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        elseif B == 12
            JwVec = plots.dualRotation.B12.JKxVec;
            KxMat = plots.dualRotation.B12.KxMat;
            for i = 1:length(KxMat(1,:))
                KxVec(i) = interp1(JwVec,KxMat(:,i), Jw, 'spline');
            end
        end
    end
end




