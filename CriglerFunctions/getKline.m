function [Kline] = getKline(nP, interpMode, B, Jw)

crieglerPlots

if strcmp(interpMode,'linear')
    if strcmp(nP,'singleRotation')
        if B == 2
            Kline = interp1(plots.singleRotation.B2.Cp,plots.singleRotation.B2.Kline,Jw,'linear', 'extrap');
        elseif B == 3
            Kline = interp1(plots.singleRotation.B3.Cp,plots.singleRotation.B3.Kline,Jw,'linear', 'extrap');
        elseif B == 4
            Kline = interp1(plots.singleRotation.B4.Cp,plots.singleRotation.B4.Kline,Jw,'linear', 'extrap');           
        end
    elseif strcmp(nP,'dualRotation')
        error('Dual rotation propeller slipstream contraction factor not available')
    end
    
end


if strcmp(interpMode,'spline')
    if strcmp(nP,'singleRotation')
        if B == 2
            Kline = interp1(plots.singleRotation.B2.Cp,plots.singleRotation.B2.Kline,Jw,'spline');
        elseif B == 3
            Kline = interp1(plots.singleRotation.B3.Cp,plots.singleRotation.B3.Kline,Jw,'spline');
        elseif B == 4
            Kline = interp1(plots.singleRotation.B4.Cp,plots.singleRotation.B4.Kline,Jw,'spline');           
        end
    elseif strcmp(nP,'dualRotation')
        error('Dual rotation propeller slipstream contraction factor not available')
    end
end

