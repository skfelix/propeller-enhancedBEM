function analysis = f_criglerPerf(inputs,geom,polares,range,delta_betadeg)

addpath('C:\Users\Felix\Documents\Mestrado\Tese\Códigos\Prop\CriglerFunctions')

Wdot = inputs.Wdot;
V = range;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
unit_system = inputs.unit_system;
x_hub = geom.x_hub;

R = D/2;

x_ref = geom.x;
b_ref  = geom.c;
beta_ref = geom.beta;

[row,~] = size(x_ref);
if row > 1
    x_ref = x_ref.';
    b_ref = b_ref.';
    beta_ref = beta_ref.';
end

% Interp
x = [0.2:0.05:0.95, 0.975, 0.99];
b = spline(x_ref,b_ref,x);
beta = spline(x_ref,beta_ref,x);
betadeg = rad2deg(beta);


delta_beta = deg2rad(delta_betadeg);
beta = beta + delta_beta;
r = x*R;
sigma = B*b./(2*pi*r);

[rho,T,patm,asound, mu] = ISA(h,unit_system);
rpm = 60*n;
omega = n*2*pi;
% J = V/(n*D);
% lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix
nP = 'singleRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline


%% Polares
meanV = min(range);
W0 = sqrt(meanV.^2 + (omega.*x.*R).^2);
Re = (rho*W0.*b/mu);
Re75 = interp1(x,Re,0.75);
meanMach = mean(W0/asound);
if meanMach > 0.45
    meanMach = 0.45;
end
if (polares == 1)
    [pol,foil] = xfoil('clark-y.dat',-10:0.5:20,Re75,0,'oper iter 150');
    pol = spera360(pol,foil);
    save('polaresNeelBright.mat','pol');
else
    load('polar-clark-y_extended.mat','pol');
end

%%

Cp = zeros(1,length(V));
Ct = zeros(1,length(V));
eta = zeros(1,length(V));
CL = zeros(1,length(x));
CLnew = zeros(1,length(x));
CD = zeros(1,length(x));

for i = 1:length(V)
    J(i) = V(i)/(n*D);
    
    W = sqrt(V(i)^2 + (omega.*x.*R).^2);
    Re = rho*W.*b/mu;
    Mach = W/340;
    
    %% Induced Velocity
    clear wbar
    Pct = Wdot/(0.5*rho*V(i)^3*pi*(D/2)^2);
    Pc = @(wbar) 2*getKappa1(nP, interpMode, B, J(i)*(1+wbar))*wbar*(1+wbar)*(1+getKappa2(nP, interpMode, B, J(i)*(1+wbar))*wbar);
    wbar = fsolve(@(wbar)Pc(wbar)-Pct,0.5);
    Pc = Pc(wbar);
    Jw(i) = J(i)*(1+wbar);
    wout(i) = wbar;
    
    k = getKappa1(nP, interpMode, B, J(i)*(1+wbar));
    e_k = getKappa2(nP, interpMode, B, J(i)*(1+wbar));
    e = e_k*k;
    
    phi = atan(1/pi*J(i)*(1+0.5*wbar)./x);
    phi0 = atan(J(i)./(pi.*x));
    x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    [KxVec1] = getKx(nP, interpMode, B, Jw(i)); % KxVec comes for x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    KxVec = interp1(x1,KxVec1,x,'linear', 'extrap');
    
    erro = 10e6;
    
    while erro > 1e-4
        
        alpha = beta - phi;
        alphadeg = rad2deg(alpha);

        for j = 1:length(x)
            CLnew(j) = spline(pol.alpha,pol.CL,alphadeg(j));
            CD(j)    = spline(pol.alpha,pol.CD,alphadeg(j));
        end
        
        erro = abs(CLnew - CL);   
        CL = CLnew;
    end

    sigmaCL = sigma.*CL; 
    eps = atan(CD./CL);
    dCqdx = pi/8*J(i)^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(sin(phi) + eps.*cos(phi)).*x.^2;
    Cq = trapz(x,dCqdx);
    dCtdx = pi/4*J(i)^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(cos(phi) - eps.*sin(phi)).*x;
    
    Ct(i) = trapz(x,dCtdx);
    Cp(i) = 2*pi*Cq;
    
    
    
    eta(i) = J(i)*Ct(i)/Cp(i);
    T(i) = Ct(i)*(rho*n^2*D^4)
    
    %%
    
%     dCqdx0 = sigmaCL.*x.^4*pi^3/8.*cos(phi).^2.*(sin(phi) + e.*cos(phi));
%     Cq0 = trapz(x,dCqdx0)
%     dCtdx0 = sigmaCL.*x.^3*pi^3/4.*cos(phi).^2.*(cos(phi) - e.*sin(phi));
%     Ct0 = trapz(x,dCtdx0)
%     Cp0 = 2*pi*Cq0;
%     
%     eta0 = J*Ct0/Cp0;
%%
    if i == 1
        [Kline] = getKline(nP, interpMode, B, Cp(i));
        FM(1) = (e +0.5*k)^(3/2)/(2*e);
        FM(2) = 0.798*Ct(i).^(3/2)./Cp(i); %sqrt(2/pi)
        FM(3) = 0.565*Ct(i).^(3/2)./(Cp(i).*Kline);
        FM(4) = 1/sqrt(2)*Ct(i).^(3/2)./Cq(i);
    end
end

% [Kline] = getKline(nP, interpMode, B, Cp);
% FM = 0.565*Ct0.^(3/2)./(Cp0*Kline)
% FM2 = (e_k*k +0.5*k)^(3/2)/(2*e_k*k)


analysis.V = V;
analysis.J = J;
analysis.Ct = Ct;
analysis.Cp = Cp;
analysis.eta = eta;
analysis.T = T;
analysis.FM = FM;
analysis.unit_system = unit_system;
