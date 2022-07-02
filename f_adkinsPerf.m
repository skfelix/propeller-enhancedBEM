function analysis = f_adkinsPerf(inputs,geom,polares,range,beta_ref_deg,x_beta_ref)

n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
unit_system = inputs.unit_system;
x_hub = geom.x_hub;

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
% c = interp1(x_ref,b_ref,x,'linear','extrap');
% beta = interp1(x_ref,beta_ref,x,'linear','extrap');
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

%
% x = geom.xi;
% c = geom.c;
% r = geom.r;
% betadeg = geom.beta_deg;
% beta = deg2rad(betadeg);
% x_hub = geom.r_hub/R;
%
% r_original = linspace(r(1),R,length(x));
% r = r_original;
% xi_original = r_original./R;     % Nondimensional radius
% c = geom.c;
% betadeg = geom.beta_deg;
%
% r(end) = interp1(xi_original,r,0.99);
% c(end) = interp1(xi_original,c,0.99);
% betadeg(end) = interp1(xi_original,betadeg,0.99);
% x = r/R;
%%
R = D/2;
r = x*R;
sigma = B*c./(2*pi*r);
r_hub = x_hub*R;

Toffset = 0;
[rho,Temp,~,asound, mu] = ISA(h,Toffset,unit_system);
rpm = 60*n;
omega = n*2*pi;
% J = V/(n*D);
% lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix

%% Polares
meanV = min(range);
W0 = sqrt(meanV.^2 + (omega.*x.*R).^2);
Re = (rho*W0.*c/mu);
Re75 = interp1(x,Re,0.75);
meanMach = mean(W0/asound);
% meanMach = 0;
if meanMach > 0.4
    meanMach = 0.4;
end

if (polares == 1)
    [pol,foil] = xfoil(geom.airfoil,-10:0.5:20,Re75,meanMach,'oper iter 150','ppar n 160','ppar t 1');
    %     pol = spera360(pol,foil);
    save('polaresNeelBright.mat','pol','foil');
end


%% main loop

V = range;

J = zeros(1,length(V));
Cp = zeros(1,length(V));
Cq = zeros(1,length(V));
Ct = zeros(1,length(V));
eta = zeros(1,length(V));
CL = zeros(1,length(x));
CD = zeros(1,length(x));

for i = length(V):-1:1

    J(i) = V(i)/(n*D);
    omega = n*2*pi;
    lambda = V(i)/(omega*R);
    chi = omega.*r/V(i);
    
    phi = atan(V(i)./(omega*r));
    
    W0 = sqrt(V(i)^2 + (omega.*x.*R).^2);
    Re = rho*W0.*c/mu;
    Mach = W0/asound;
    
    W = W0; % Initial Guess
    
    erro = 10e6;
    while erro > 1e-5
        
        phit = atan(tan(phi(end)).*x(end));
        f  = 0.5*B*(1-x)./sin(phit); % adkins
        F  = 2/pi*acos(exp(-f));
        fh = 0.5*B*(r-r_hub)./(r_hub*sin(phi));
        Fh = 2/pi*acos(exp(-fh));
        F  = F.*Fh;
        
        alpha = beta - phi;
        alphadeg = rad2deg(alpha);
        
        pol360 = spera360(pol,foil,geom,1);
        for j = 1:length(r)
            %             pol_3D = postStallDuSelig(pol,R,V(i),n,r(j),c(j));
            %             pol3D = postStallCorriganSchillings(pol,r(j),c(j));
%             foil.max_thickness = 0.12;
%             pol360 = spera360(pol,foil,geom,1);
            %             pol360 = viternaCorrigan360(pol3D,foil);
            %             plot(pol.alpha,pol.CL,pol3D.alpha,pol3D.CL,pol360.alpha,pol360.CL)
            CL(j) = interp1(pol360.alpha,pol360.CL,alphadeg(j),'linear','extrap');
            CD(j) = interp1(pol360.alpha,pol360.CD,alphadeg(j),'linear','extrap');
        end
        % CD Mach Reynolds Correction, using Re75 as reference. Hernandez &
        % Crespo (1987)
        Re = rho*W.*c/mu;
        Mach = W/asound;
        CD = CD.*(Re75./Re).^0.11; % Relation from experimental data in NACA TR 586. More reliable because BL is not always turbulent on airfoil
        CD = CD.*sqrt(1-meanMach.^2)./sqrt(1-Mach.^2);
        
        Cy = CL.*cos(phi) - CD.*sin(phi);
        Cx = CL.*sin(phi) + CD.*cos(phi);
        
        Ka = Cy./(4*sin(phi).^2);
        Kb = Cx./(4*cos(phi).*sin(phi));
        
        a = sigma.*Ka./(F - sigma.*Ka);
        b = sigma.*Kb./(F + sigma.*Kb);
        
        if phi(end) < 0
            phi(end) = deg2rad(1);
        end
        
        phinew = atan(V(i).*(1+a)./(omega.*r.*(1-b)));
        
        erro = max(abs(phi-phinew));
        k = 0.4;
        phi = (1-k)*phi + k*phinew;
    end
    
    W = V(i).*(1+a)./sin(phi);
    
    dT = 0.5*rho.*W.^2*B.*c.*Cy;
    T(i) = trapz(r,dT)*0.95 %%%%%%%%%%%%%%%
    dQ = 0.5*rho.*W.^2*B.*c.*Cx.*r;
    Q(i) = trapz(r,dQ);
    P(i) = omega*Q(i); % 2*pi*n*Q
    Pel(i) = P(i)/0.85; 
    Ct(i) = T(i)./(rho*n^2*D^4);
    Cq(i) = Q(i)./(rho*n^2*D^5);
    Cp(i) = P(i)./(rho*n^3*D^5);
    eta(i) = J(i)*Ct(i)/Cp(i);
end

analysis.rpm = rpm;
analysis.n = n;
analysis.V = V(Ct~=0);
analysis.J = J(Ct~=0);
analysis.Ct = Ct();
analysis.Cp = Cp(Ct~=0);
analysis.Cq = Cq(Ct~=0);
analysis.eta = eta(Ct~=0);
analysis.T = T(Ct~=0);
analysis.P = P(Ct~=0);
analysis.Pel = Pel(Ct~=0);
analysis.Q = Q(Ct~=0);
analysis.unit_system = unit_system;
