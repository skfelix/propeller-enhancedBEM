function plotGeom(geom,title,lgd)

% fontBold = 'bold';
fontSize = 12;
lgdSize = 10;
solidLine = '-';
lineWidth = 1;
% lineVec = {'-','--',':','-.'};
colorVec = {'k','r','b','g','m','c','y'};
markerVec = {'o','s','^','d','v','+','x'};
% markerNumber = 10;

lgd_pitch = lgd;

sgtitle(title);
%% Chord
subplot(1,2,1);
% title('\textbf{a)}','Interpreter','Latex','FontSize',fontSize);
hold on
for i = 1:length(geom)
    if strcmp(geom{1}.unit_system,'imperial')
        plot(geom{i}.x,geom{i}.c*12,...
        'LineStyle',solidLine,...
        'Color',colorVec{i},...
        'Marker',markerVec{i},...
        'LineWidth',lineWidth)
    else
        plot(geom{i}.x,geom{i}.c*1e3,...
        'LineStyle',solidLine,...
        'Color',colorVec{i},...
        'Marker',markerVec{i},...
        'LineWidth',lineWidth)
    end
end
xlabel('x = r/R','Interpreter','Latex','FontSize',fontSize);
if strcmp(geom{1}.unit_system,'imperial')
    ylabel('c [ft]','Interpreter','Latex','FontSize',fontSize);
else
    ylabel('c [mm]','Interpreter','Latex','FontSize',fontSize);
end
grid on; grid minor;
if ~isempty(lgd)
    legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
end

%% Beta
subplot(1,2,2);
% title('\textbf{b)}','Interpreter','Latex','FontSize',fontSize);
hold on
for i = 1:length(geom)
    plot(geom{i}.x,geom{i}.betadeg,...
        'LineStyle',solidLine,...
        'Color',colorVec{i},...
        'Marker',markerVec{i},...
        'LineWidth',lineWidth)
    lgd_pitch{i} = ['p75 = ' num2str(geom{i}.pitch75inch,'%.1f') '"'];
end
xlabel('x = r/R','Interpreter','Latex','FontSize',fontSize);
ylabel('$\beta$ [deg]','Interpreter','Latex','FontSize',14);
grid on; grid minor;
legend(lgd_pitch,'Interpreter','Latex','Location','best','FontSize',lgdSize) 


% ylim([0 inf])