function plotPerf(perf,title,lgd,Cp1,varargin)

if nargin < 5
    plotVelocity = 0;
else 
    plotVelocity = varargin{1};
end

if nargin < 6
    plotThrust = 0;
else 
    plotThrust = varargin{2};
end

if nargin < 7
    plotStatic = 0;
else 
    plotStatic = varargin{3};
end


% fontBold = 'bold';
fontSize = 12;
lgdSize = 12;
solidLine = '-';
lineWidth = 1;
% lineVec = {'-','--',':','-.'};
colorVec = {'k','r','b','g','m','c','#D95319','#7E2F8E'};
markerColorVec = {'k','r','b','g','m','c','#D95319','#7E2F8E'};
markerVec = {'o','s','^','d','v','+','x','*'};
markerQuantity = 10;

maxCpVec = zeros(1,length(perf));
maxCtVec = zeros(1,length(perf));
maxTVec = zeros(1,length(perf));
xMaxVec = zeros(1,length(perf));

if plotVelocity
    xVec = 'V';
else
    xVec = 'J';
end

for i = 1:length(perf)
    xMaxVec(i) = max(perf{i}.(xVec));
end
xMax = max(xMaxVec);

if plotThrust ==0
    sgtitle(title)
    %% Efficiency
    subplot(2,2,1);
    hold on
    for i = 1:length(perf)
        plot(perf{i}.(xVec),perf{i}.eta,...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
    end
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    ylabel('$\eta$','Interpreter','Latex','FontSize',fontSize);
    grid on; grid minor;
    % legend(lgd,'Interpreter','Latex','Location','best')
    xLim = xlim;
    xlim([0 xLim(2)]);
    ylim([0 1]);
    
    %% Power coefficient - Cp
    subplot(2,2,2);
    hold on
    for i = 1:length(perf)
        plot(perf{i}.(xVec),perf{i}.Cp/0.85,...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxCpVec(i) = max(perf{i}.Cp);
    end
    if Cp1 ~= 0
        plot([0 xMax],[Cp1 Cp1],'k--','lineWidth',2)
        maxCp = max([maxCpVec Cp1])*1.1;
    else
        maxCp = max(maxCpVec)*1.1;
    end
    
    % maxCp = 1.1*Cp1;
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    ylabel('$C_P$','Interpreter','Latex','FontSize',14);
    grid on; grid minor;
    % legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
    axis([0 xMax 0 maxCp])
    
    %% Thrust Coefficient - Ct
    subplot(2,2,3);
    hold on
    j = 1;
    for i = 1:length(perf)
        plot(perf{i}.(xVec),perf{i}.Ct,...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxCtVec(i) = max(perf{i}.Ct);
%         if sum(perf{i}.Cp > 0.95*Cp1)>0 .... trocar i por j = j + 1
%                 V_Cp_max(j) = interp1(perf{i}.Cp,perf{i}.V,Cp1,'linear','extrap');
%             Ct_Cp_max(j) = interp1(perf{i}.(xVec),perf{i}.Ct,V_Cp_max(j),'linear','extrap');
%             j = j + 1;
%         end
    end
    if j > 1; plot(V_Cp_max,Ct_Cp_max,'--k','lineWidth',2); end
    maxCt = max(maxCtVec)+0.02;
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    ylabel('$C_T$','Interpreter','Latex','FontSize',14);
    grid on; grid minor;
    % legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
    axis([0 xMax 0 maxCt]);
    
    %% Thrust force - T
    subplot(2,2,4);
    % title('\textbf{b)}','Interpreter','Latex','FontSize',fontSize);
    hold on
    j = 1;
    for i = 1:length(perf)
        plot(perf{i}.(xVec),perf{i}.T,...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxTVec(i) = max(perf{i}.T);
%         if sum(perf{i}.Cp > 0.95*Cp1)>0 .... trocar i por j = j + 1
%                 V_Cp_max(j) = interp1(perf{i}.Cp,perf{i}.V,Cp1,'linear','extrap');
%             T_Cp_max(j) = interp1(perf{i}.(xVec),perf{i}.T,V_Cp_max(j),'linear','extrap');
%             j = j + 1;
%         end
    end
    if j > 1; plot(V_Cp_max,T_Cp_max,'--k','lineWidth',2); end
    maxT = max(maxTVec)*1.25;
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    if strcmp(perf{end}.unit_system,'imperial')
        ylabel('$T$ [lbf]','Interpreter','Latex','FontSize',14);
    else
        ylabel('$T$ [N]','Interpreter','Latex','FontSize',14);
    end
    grid on; grid minor;
    lgd{end+1} = '$C_{p_{max}}$';
    legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize,'AutoUpdate','off')
    axis([0 xMax 0 maxT]);
    
elseif plotThrust == 1  
    
    sgtitle(title)
    subplot(121)
    hold on
    j = 1;
    for i = 1:length(perf)
        plot(perf{i}.(xVec),movmean(perf{i}.T/9.81,4), ...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxTVec(i) = max(perf{i}.T/9.81);
    end
    maxT = max(maxTVec)*1.25;
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    elseif plotStatic
        xlabel('$RPM$','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    if strcmp(perf{end}.unit_system,'imperial')
        ylabel('$T$ [lbf]','Interpreter','Latex','FontSize',14);
    else
        ylabel('$T$ [kgf]','Interpreter','Latex','FontSize',14);
    end
    grid on; grid minor;
    lgd{end+1} = '$C_{p_{max}}$';
    legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize,'AutoUpdate','off')
    axis([0 xMax 0 maxT]);
    
    subplot(122)
    hold on
    for i = 1:length(perf)
        plot(perf{i}.(xVec),movmean(perf{i}.P/0.85,4),...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxPVec(i) = max(perf{i}.P);
    end
    
    maxP = max(maxPVec)*1.5;
    
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    elseif plotStatic
        xlabel('$RPM$','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    ylabel('$P$ [W]','Interpreter','Latex','FontSize',14);
    grid on; grid minor;
    % legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
    axis([0 xMax 0 maxP])
    
    %%%%%%%%%%%
elseif plotThrust ==2
    sgtitle(title)
    subplot(121)
    hold on
    j = 1;
    for i = 1:length(perf)
        plot(perf{i}.(xVec),movmean(perf{i}.T,4), ...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxTVec(i) = max(perf{i}.T);
    end
    maxT = max(maxTVec)*1.25;
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    elseif plotStatic
        xlabel('$RPM$','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    if strcmp(perf{end}.unit_system,'imperial')
        ylabel('$T$ [lbf]','Interpreter','Latex','FontSize',14);
    else
        ylabel('$T$ [kgf]','Interpreter','Latex','FontSize',14);
    end
    grid on; grid minor;
    lgd{end+1} = '$C_{p_{max}}$';
    legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize,'AutoUpdate','off')
    axis([0 xMax 0 maxT]);
    
    subplot(122)
    hold on
    for i = 1:length(perf)
        plot(perf{i}.(xVec),movmean(perf{i}.Q,4),...
            'LineStyle',solidLine, ...
            'Color',colorVec{i}, ...
            'Marker',markerVec{i}, ...
            'MarkerEdgeColor', markerColorVec{i}, ...
            'LineWidth',lineWidth, ...
            'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
        maxQVec(i) = max(perf{i}.Q);
    end
    
    maxQ = max(maxQVec)*1.5;
    
    if plotVelocity
        xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
    elseif plotStatic
        xlabel('$RPM$','Interpreter','Latex','FontSize',fontSize);
    else
        xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
    end
    ylabel('$Q$ [N.m]','Interpreter','Latex','FontSize',14);
    grid on; grid minor;
    % legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
    axis([0 xMax 0 maxQ])
end