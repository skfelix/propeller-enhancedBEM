function analysis = f_criglerPerfMod(inputs,geom,polares,range,beta_ref_deg,x_beta_ref,varargin)

if nargin < 7
    regime = 'dynamic';
else 
    regime = varargin{1};
end

addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\CriglerFunctions')

V = range;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
unit_system = inputs.unit_system;
x_hub = geom.x_hub;

R = D/2;

% Original geometry
x_original = geom.x;
c_original  = geom.c;
beta_original = geom.beta;
betadeg_original = geom.betadeg;

[row,~] = size(x_original);
if row > 1
    x_original          = x_original.';
    c_original          = c_original.';
    beta_original       = beta_original.';
    betadeg_original    = betadeg_original.';
end

% Interpolating to default x vector
x = [0.2:0.05:0.95, 0.975, 0.99];
c = spline(x_original,c_original,x);
beta = spline(x_original,beta_original,x);
betadeg = spline(x_original,betadeg_original,x);

% Matching beta ref for r = 0.75R or 0.70R
delta_betadeg = beta_ref_deg - interp1(x,betadeg,x_beta_ref);
delta_beta = deg2rad(delta_betadeg);
beta_matching = beta + delta_beta;
betadeg_matching = betadeg + delta_betadeg;

beta = beta_matching;
betadeg = betadeg_matching;

r = x*R;
sigma = B*c./(2*pi*r);

T_offset = 0;
[rho,~,~,asound, mu] = ISA(h,T_offset,unit_system);
rpm = 60*n;
omega = n*2*pi;
S = pi*R^2; % Projected area of the helix
nP = 'singleRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'linear'; % 'linear' or 'spline

AF = 1e5/(32*R^5)*trapz(x*R, c.*(x*R).^3);
TAF = B*AF;
sigma75 = B*interp1(x,c,0.75)/(2*pi*0.75*R);

%% Polares
meanV = min(range);
W0 = sqrt(meanV.^2 + (omega.*x.*R).^2);
Re = (rho*W0.*c/mu);
Re75 = interp1(x,Re,0.75);
meanMach = mean(W0/asound);
if meanMach > 0.45
    meanMach = 0.45;
end
if (polares == 1)
    [pol,foil] = xfoil(geom.airfoil,-4:0.5:19.5,Re75,meanMach,'oper iter 150','ppar n 200','ppar t 1');
%     pol = spera360(pol,foil);
    save('polaresNeelBright.mat','pol','foil');
elseif polares == 0
    load('polaresNeelBright.mat','pol','foil');
elseif polares == 2
    load('polaresBroeren_interp.mat','polaresBroeren_interp');
    pol = polaresBroeren_interp.clean;
elseif polares == 3
    load('polaresBroeren_interp.mat','polaresBroeren_interp');
    pol = polaresBroeren_interp.roughness1;
elseif polares == 4
    load('polaresBroeren_interp.mat','polaresBroeren_interp');
    pol = polaresBroeren_interp.streamwise1;
elseif polares == 5
    load('polaresBroeren_interp.mat','polaresBroeren_interp');
    pol = polaresBroeren_interp.horn;
elseif polares == 6
    load('polaresBroeren_interp.mat','polaresBroeren_interp');
    pol = polaresBroeren_interp.spanwiseRidge;
end

%%

Cp = zeros(1,length(V));
Ct = zeros(1,length(V));
eta = zeros(1,length(V));
Jw = zeros(1,length(V));
CL = zeros(1,length(x));
CD = zeros(1,length(x));

pol360 = spera360(pol,foil,geom,1);

wbar = 0.1;
breakOutFor = 0;
% for i = 1:length(V)
for i = length(V):-1:1
    
    if breakOutFor
        break;
    end
    
    breakOutFor = 0;
    
    J(i) = V(i)/(n*D);
    
%     phi0 = atan(V(i)./(omega*r));
    phi = atan(J(i)*(1+0.5*wbar)./(pi*x));    
    W0 = sqrt(V(i)^2 + (omega.*x.*R).^2);
    
    erro = 10e6;
    count = 1;
    while erro > 1e-2
        
        alpha = beta - phi;
        alphadeg = rad2deg(alpha);
        
        if sum(alphadeg > max(pol.alpha)) > round(length(r)/2) && polares > 1
%             breakOutFor = 1;
        end
        
%         pol360 = spera360(pol,foil,geom,1);
        for j = 1:length(r)
% %             pol3D = postStallDuSelig(pol,R,V(i),n,r(j),c(j));
%             pol3D = postStallCorriganSchillings(pol,r(j),c(j));
%             foil.max_thickness = 0.12;            
%             pol360 = spera360(pol,foil,geom,1);
%             pol360 = pol;
%             pol360 = viternaCorrigan360(pol3D,foil);
%             plot(pol.alpha,pol.CL,pol360.alpha,pol360.CL)
            CL(j) = interp1(pol360.alpha,pol360.CL,alphadeg(j),'linear','extrap');
            CD(j) = interp1(pol360.alpha,pol360.CD,alphadeg(j),'linear','extrap');
        end
        
        % CD Mach Reynolds Correction, using Re75 as reference. Hernandez &
        % Crespo (1987)
        W = V(i)*(1+0.5*wbar*cos(phi).^2)./sin(phi);
        Re = rho*W.*c/mu;
        Mach = W/asound;   
%         CD = CD.*(Re75./Re).^0.20; % Relation form turbulent BL eq: Cf = 0.074/Re^(1/5)
        CD = CD.*(Re75./Re).^0.11; % Relation from experimental data in NACA TR 586. More reliable because BL is not always turbulent on airfoil
        CD = CD.*sqrt(1-meanMach.^2)./sqrt(1-Mach.^2);
        
        Cy = CL.*cos(phi) - CD.*sin(phi);
        Cx = CL.*sin(phi) + CD.*cos(phi);

        dCqdx = 1/8*pi*x.^2*J(i)^2.*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigma.*Cx;
        if strcmp(regime,'static')
            Cq(i) = trapz(x,dCqdx)*1.25;
        else
            Cq(i) = trapz(x,dCqdx);
        end
        
        Cp(i) = 2*pi*Cq(i);
        
        Pct = 64/pi*(n*R/V(i))^3*Cp(i); % Cp to Pc
        options = optimset('Display','off');
        Pc = @(wbar) 2*getKappa1(nP, interpMode, B, J(i)*(1+wbar))*wbar*(1+wbar)*(1+getKappa2(nP, interpMode, B, J(i)*(1+wbar))*wbar);
        wbarnew = fsolve(@(wbar)Pc(wbar)-Pct,wbar,options);
        Jw(i) = J(i)*(1+wbarnew);
        
        phinew = atan(J(i)*(1+0.5*wbarnew)./(pi*x));   
        
        erro = max(abs(wbar-wbarnew));
%         erro = max(abs(phi-phinew));
        
        wbar = wbarnew;
        phi = phinew;
        
        count = count + 1;
        if count > 50
            break;
        end
    end
    W0_out{i} = W0;
    alpha_out{i} = alphadeg;
    Re_out{i} = Re;
    Re75_out{i} = Re75;
    Mach_out{i} = Mach;
    meanMach_out{i} = meanMach;
    
    k = getKappa1(nP, interpMode, B, Jw(i));
    e_k = getKappa2(nP, interpMode, B, Jw(i));
    e = e_k*k;
    dCtdx = 1/4*pi*x*J(i)^2.*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigma.*Cy;
    if strcmp(regime,'static')
            Ct(i) = trapz(x,dCtdx)*0.95;
        else
            Ct(i) = trapz(x,dCtdx)*0.85;
    end
    eta(i) = J(i)*Ct(i)/Cp(i);
    T(i) = Ct(i)*(rho*n^2*D^4)
    Q(i) = Cq(i)*(rho*n^2*D^5);
    P(i) = Cp(i)*(rho*n^3*D^5);
    
    %%
%     if i == 1 
%         [Kline] = getKline(nP, interpMode, B, Cp(i));
%         FM(1) = (e +0.5*k)^(3/2)/(2*e);
%         FM(2) = 0.798*Ct(i).^(3/2)./Cp(i); %sqrt(2/pi)
%         FM(3) = 0.565*Ct(i).^(3/2)./(Cp(i).*Kline);
%         analysis.FM = FM;
%     end
end

analysis.rpm = rpm;
analysis.n = n;
analysis.V = V(Ct~=0);
analysis.J = J(Ct~=0);
analysis.Ct = Ct(Ct~=0);
analysis.Cp = Cp(Ct~=0);
analysis.Cpel = analysis.Cp/0.85;
analysis.Cq = Cq(Ct~=0);
analysis.eta = eta(Ct~=0);
analysis.T = T(Ct~=0);
analysis.P = P(Ct~=0);
analysis.Pel = analysis.P/0.85;
analysis.Q = Q(Ct~=0);
analysis.unit_system = unit_system;
analysis.Re75 = Re75;
analysis.Mach75 = meanMach;

% For Cp slice
% analysis.x = x;
% analysis.alpha = alpha_out;
% analysis.W0 = W0_out;
% analysis.Re = Re_out;
% analysis.Re75 = Re75_out;
% analysis.Mach = Mach_out;
% analysis.meanMach = meanMach_out;