% %% MAIN PROP - ALL BETA
% clc; fclose all;
% addpath('C:\Users\Felix\Documents\ProgramasAero\xfoil');
% addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop')
% addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\AuxiliarFunctions');
% addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\CriglerFunctions');
% addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\GeometryAndPerfData');
% %%%% Baseline Mestrado - metric
% Wdot = 500*0.65;
% Wdot_max = 500;
% % V = 15; % set1
% V = 20; %25 ok % 20 too low
% D = 0.35;
% B = 2;
% h = 1100;
% unit_system = 'metric';
% airfoil = 'NACA23012';
% x_hub = 0.10;
% 
% inputs.Wdot = Wdot; inputs.Wdot_max = Wdot_max; inputs.V = V; inputs.D = D; inputs.B = B;
% inputs.h = h; inputs.unit_system = unit_system; inputs.airfoil = airfoil; inputs.x_hub = x_hub;
% 
% rho = ISA(h,0,unit_system);
% polares  = 1;
% propType = 3;
% x        = [0.2:0.05:0.95, 0.975, 0.99];
% 
% n = [8 9 10 11 12]*1e3/60;
% 
% roda = 0;
% if roda
%     % Design
%     for i = 1:length(n)
%         inputs.n = n(i);
%         desRot(i) = f_adkinsDesign(inputs,polares,propType,0.5,x);
%         T_des(i) = desRot(i).T;
%     end
%     
%     % Performance
%     x_beta_ref = 0.75;
%     range = 1:2:25;
%     inputs.n = 4000/60;
%     for i = 1:length(n)
%         beta_ref_deg = interp1(desRot(i).x,desRot(i).betadeg,x_beta_ref);
%         anRot(i) = f_criglerPerfMod(inputs,desRot(i),polares,range,beta_ref_deg,x_beta_ref);
%     end
%     save('mainV3_nSweep.mat','desRot','anRot')
% else
%     load('mainV3_nSweep.mat')
% end
% 
% 
% %% Plots
% clear geom perf
% nfig = 1;
% rho = ISA(h);
% for i = 1:length(n)
%     lgd{i} = string(n(i)*60) + ' ' + 'rpm';
%     geom{i} = desRot(i);
%     perf{i} = anRot(i);
%     Cp1(i) = Wdot_max/(rho*n(i)^3*D^5);
% end
% Cp1 = max(Cp1);
% title = strcat('Design RPM variation');
% set(groot,'defaultAxesTickLabelInterpreter','latex'); 
% 
% % Geometry
% close all
% figure(nfig);
% plotGeom(geom,title,lgd)
% 
% % Perf
% nfig = nfig+1;
% plotVelocity = 1;
% plotThrust = 1;
% figure(nfig);
% plotPerf(perf,title,lgd,Cp1,plotVelocity,plotThrust)

%% MAIN PROP - CL sweep
clc; fclose all; 
addpath('C:\Users\Felix\Documents\ProgramasAero\xfoil');
addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop')
addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\AuxiliarFunctions');
addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\CriglerFunctions');
addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\GeometryAndPerfData');
%%%% Baseline Mestrado - metric
Wdot_design = 500*0.55;
Wdot_max = 500;
V = 18;
n = 8e3/60;
D = 0.33;
B = 2;
h = 1100;
unit_system = 'metric';
airfoil = 'NACA23012';
x_hub = 0.10;

inputs.Wdot_max = Wdot_max; inputs.Wdot = Wdot_design; inputs.V = V; inputs.n = n; inputs.D = D; inputs.B = B; 
inputs.h = h; inputs.unit_system = unit_system; inputs.airfoil = airfoil; inputs.x_hub = x_hub;

rho = ISA(h,0,unit_system);
polares  = 1;
propType = 3;
x        = [0.15:0.05:0.95, 0.975, 0.99];

Cp1 = Wdot_max/(rho*n^3*D^5);

CL_design = [0.50 0.515 0.53 0.56];

roda = 1;
if roda
    % Design
    for i = 1:length(CL_design)
        desCL(i) = f_adkinsDesign(inputs,polares,propType,CL_design(i),x);
        T_des(i) = desCL(i).T;
    end
%     desCL(6) = [];
    desCL(5) = desCL(2);
    vec = [0.90 0.82 0.82 0.87 0.92 0.96 0.98 1 1 1 1 1 1 1 1 1 1 1 1];
    desCL(5).c = desCL(5).c.*vec;
    CL_design(5) = CL_design(2);
    
    % Performance
    x_beta_ref = 0.75;
    range = .1:1:15;
    inputs.n = 3000/60;
    for i = 1:length(desCL)
        beta_ref_deg = interp1(desCL(i).x,desCL(i).betadeg,x_beta_ref);
        anCL(i) = f_criglerPerfMod(inputs,desCL(i),polares,range,beta_ref_deg,x_beta_ref);
    end
    save('mainV3_CLSweep.mat','desCL','anCL')
else
    load('mainV3_CLSweep.mat')
end

%%
% MSc_v3.geom = desCL(5);
% MSc_v3.perf = anCL(5);
% save('MSc_v3.mat','MSc_v3')

%% Plots
clear geom perf
nfig = 2;
for i = 1:length(anCL)
    lgd{i} = ['$C_L = ' num2str(CL_design(i),'%.2f') '$'];
        geom{i} = desCL(i);
    perf{i} = anCL(i);
end
% lgd{5} = ['$C_L = ' num2str(0.52,'%.2f') '$ Mod. Hub'];
% geom{5} = desCL(5);
% perf{5} = anCL(5);
rho = ISA(h);
Cp1 = Wdot_max/(rho*n^3*D^5);

title = strcat('Design CL variation');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

% Geometry
close all
figure(nfig);
plotGeom(geom,title,lgd)

%% Perf
nfig = nfig+1;
plotVelocity = 1;
figure(nfig);
plotPerf(perf,title,lgd,Cp1,plotVelocity,1)

%% The chosen propeller anCL(5)

Wdot_design = 500*0.55;
Wdot_max = 500;
V = 18;
n = 8e3/60;
D = 0.33;
B = 2;
h = 1100;
unit_system = 'metric';
airfoil = 'NACA23012';
x_hub = 0.10;

inputs.Wdot_max = Wdot_max; inputs.Wdot = Wdot_design; inputs.V = V; inputs.n = n; inputs.D = D; inputs.B = B; 
inputs.h = h; inputs.unit_system = unit_system; inputs.airfoil = airfoil; inputs.x_hub = x_hub;

rho = ISA(h,0,unit_system);
polares  = 1;
propType = 3;
x        = [0.15:0.05:0.95, 0.975, 0.99];

Cp1 = Wdot_max/(rho*n^3*D^5);

CL_design = [0.515];
rpm = [2500 3000 3500];

roda = 1;
if roda
    % Design
    desOfficial = f_adkinsDesign(inputs,polares,propType,CL_design,x);
    T_des = desOfficial.T;
    vec = [0.90 0.82 0.82 0.87 0.92 0.96 0.98 1 1 1 1 1 1 1 1 1 1 1 1];
    desOfficial.c = desOfficial.c.*vec;
    
    % Performance
    x_beta_ref = 0.75;
    range{1} = .1:0.5:10;
    range{2} = .1:0.5:12;
    range{3} = .1:0.5:14;
    for i = 1%:length(rpm)
        inputs.n = rpm(i)/60;
        beta_ref_deg = interp1(desOfficial.x,desOfficial.betadeg,x_beta_ref);
        anOfficial(i) = f_criglerPerfMod(inputs,desOfficial,polares,range{i},beta_ref_deg,x_beta_ref);
    end
    save('MScV3_official.mat','desOfficial','anOfficial')
else
    load('MScV3_official.mat')
end

%% Plots
clear geom perf
nfig = 2;
for i = 1:length(rpm)
    lgd{i} = ['RPM = ' num2str(rpm(i),'%.0f')];
    geom = {desOfficial};
    perf{i} = anOfficial(i);
end
rho = ISA(h);
Cp1 = Wdot_max/(rho*n^3*D^5);

title = strcat('Design CL variation');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

% Geometry
close all
figure(nfig);
plotGeom(geom,title,lgd)

% Perf
nfig = nfig+1;
plotVelocity = 1;
figure(nfig);
plotPerf(perf,title,lgd,Cp1,plotVelocity,2)

