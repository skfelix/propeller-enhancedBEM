function out = f_adkinsDesign(inputs,polares,propType,CL_design,x)
addpath('C:\Users\Felix\Documents\ProgramasAero\xfoil')
% Inputs
Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
airfoil = inputs.airfoil;
unit_system = inputs.unit_system;
x_hub = inputs.x_hub;

T_offset = 0;
[rho,~,~,asound, mu] = ISA(h,T_offset,unit_system);
R = D/2;
rpm = 60*n;
omega = 2*pi*n;
J = V/(n*D);
lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix

Cp_design = Wdot/(rho*n^3*D^5);
Pc = Wdot/(0.5*rho*V^3*pi*R^2);
%% Polares - Executar somente se mudar as condições de contorno
W0 = sqrt(V.^2 + (omega.*x.*R).^2);
Re = (rho*W0.*(0.075*D)/mu);
Re75 = interp1(x,Re,0.75);

%  polares = 1 ; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    %     [pol foil] = xfoil('clark-y.dat',-5:0.1:15,mean(Re),0.3,'oper iter 300');
%     [p, f] = xfoil('clark-y.dat',0:0.5:10,0.75e6,0,'oper iter 150');
    [pol,foil] = xfoil(airfoil,-4:0.5:12,Re75,0,'oper iter 150','ppar n 160','ppar t 1');
%         p = spera360(p,f);
    for i = 1:length(x)
        polRe{i} = pol;
        foilRe{i} = foil;
    end
    %     save('polares.mat','polRe','foilRe','pol','foil');
%     save('polar_clarky.mat','polRe','foilRe');
else
    load('polar_clarky.mat');
end

% %% Polares - Executar somente se mudar as condições de contorno
% % polares = 1 ; % 0 para não gerar polar, 1 para gerar
% if (polares == 1)
%     %     [pol foil] = xfoil('clark-y.dat',-5:0.1:15,mean(Re),0.3,'oper iter 300');
%     for i = 1:length(x)
%         [p f] = xfoil('clark-y.dat',-2:1:11,Re(i),Mach(i),'oper iter 150');
%         polRe{i} = p;
%         foilRe{i} = f;
%     end
%     %     save('polares.mat','polRe','foilRe','pol','foil');
%     save('polar_clarky.mat','polRe','foilRe');
% else
%     load('polar_clarky.mat');
% end

%% Prandtl Factor

r = x*R;
W = sqrt(V^2+(omega*r).^2);
Mach = W/asound;

zeta = 0.2;
erro = inf;

while erro > 1e-4 
    
    phit = atan(lambda*(1+0.5*zeta));
    
    chi = omega*r/V;
    f = 0.5*B*(1-x)/sin(phit);
    F = 2/pi*acos(exp(-f));
    phi = atan((1+0.5*zeta).*lambda./x);
    phideg = rad2deg(phi);
    Fh = 2/pi*acos(exp(-B/2*(x-x_hub)./(x_hub*sin(phi))));
    F = F.*Fh;
    G = F.*chi.*cos(phi).*sin(phi);

    
    %% Prop Type
    
    if propType == 1 % max CL/CD
        for i = 1:length(r)
            cl_cd_max = max(polRe{i}.CL_CD);
            alphadeg(i) = interp1(polRe{i}.CL_CD,polRe{i}.alpha,cl_cd_max);
            CL(i) = interp1(polRe{i}.alpha,polRe{i}.CL,alphadeg(i));
            CD(i) = interp1(polRe{i}.alpha,polRe{i}.CD,alphadeg(i));
        end
        e = 1/cl_cd_max;
    elseif propType == 2 % min CD
        for i = 1:length(r)
            cd_min = min(polRe{i}.CD);
            alphadeg(i) = interp1(polRe{i}.CD,polRe{i}.alpha,cd_min);
            CL(i) = interp1(polRe{i}.alpha,polRe{i}.CL,alphadeg(i));
            CD(i) = cd_min;
        end
        e = CD./CL;
    elseif propType == 3 % fixed CL | CL design
        for i = 1:length(r)
            alphadeg(i) = interp1(polRe{i}.CL,polRe{i}.alpha,CL_design);
            CL(i) = CL_design;
            CD(i) = interp1(polRe{i}.alpha,polRe{i}.CD,alphadeg(i));
        end
        e = CD./CL;
    elseif propType == 4 % Takeoff
        for i = 1:length(r)
            alphadeg(i) = interp1(polRe{i}.CL,polRe{i}.alpha,CL_design);
            CL(i) = CL_design;
            CD(i) = interp1(polRe{i}.alpha,polRe{i}.CD,alphadeg(i));
        end
        e = CD./CL;
    end
    
    alpha = deg2rad(alphadeg);
    Cy = CL.*(cos(phi) - e.*sin(phi));
    Cx = CL.*(sin(phi) + e.*cos(phi));
    
    %%
    
    Wc = 4*pi*lambda*G*V*R*zeta./(CL*B);
    Re = rho*Wc/mu;
    
    a = 0.5*zeta*cos(phi).^2.*(1-e.*tan(phi));
    b = 0.5*zeta./chi.*cos(phi).*sin(phi).*(1+e./tan(phi)); 
%     a(end) = 0;
%     b(end) = 0;
    W = V*(1+a)./sin(phi);
    Mach = W/asound;
    c = Wc./W;
    beta = alpha + phi;
    betadeg = rad2deg(beta);
    
    dI1 = 4*x.*G.*(1-e.*tan(phi));
    dI2 = lambda.*(0.5*dI1./x).*(1+e./tan(phi)).*sin(phi).*cos(phi);
    dJ1 = 4*x.*G.*(1+e./tan(phi));
    dJ2 = 0.5*dJ1.*(1-e.*tan(phi)).*cos(phi).^2;
    
    I1 = trapz(x,dI1);
    I2 = trapz(x,dI2);
    J1 = trapz(x,dJ1);
    J2 = trapz(x,dJ2);
    zetanew = -(0.5*J1./J2) + sqrt((0.5*J1./J2).^2 + Pc./J2);

    
    Tc = I1.*zeta - I2.*zeta.^2;
    
    erro = abs(zeta - zetanew);
    
    zeta = zetanew;
    
end

Jw = mean(V*(1+a)/n/D);
sigma = B*c./(2*pi.*r);
sigmaCL = sigma.*CL;
cCL = c.*CL;

eta = Tc/Pc;
dT = 0.5*rho.*W.^2*B.*c.*Cy;
T = trapz(r,dT);
dQ = 0.5*rho.*W.^2*B.*c.*Cx.*r;
Q = trapz(r,dQ);
P = omega*Q;
Ct = T/(rho*n^2*D^4);
Cp = P/(rho*n^3*D^5);

% Geometry

beta75 = interp1(x,beta,0.75);
betadeg75 = rad2deg(beta75);
sigma75 = B*interp1(x,c,0.75)/(2*pi*0.75*R);
AF = 1e5/(32*R^5)*trapz(x*R, c.*(x*R).^3);
TAF = B*AF;
pitch = 2*pi*(x*R).*tan(beta);
pitch75 = 2*pi*(0.75*R)*tan(beta75);

%% Outputs

out.pitch = pitch;
out.pitch75 = pitch75;
if strcmp(unit_system,'imperial')
    out.pitch75inch = pitch75*12;
else
    out.pitch75inch = pitch75/0.305*12;
end
out.x = x;
out.r = r;
out.x_hub = x_hub;
out.c = c;
out.cCL = cCL;
out.beta = beta;
out.betadeg = betadeg;
out.x_ref = 0.75;
out.beta75 = beta75;
out.betadeg75 = betadeg75;
out.phideg = rad2deg(phi);
out.phi = phi;
out.alpha = alpha;
out.alphadeg = alphadeg;
out.sigma = sigma;
out.sigmaCL = sigmaCL;
out.cCL = cCL;
out.Ct = Ct;
out.Cp = Cp;
out.CL = CL;
out.Cq = Cp/omega;
out.T = T;
out.P = P;
out.Q = Q;
out.J = J;
out.eta = eta;
out.sigma75 = sigma75;
out.TAF = TAF;
out.airfoil = airfoil;
% an = crieglerSingle_perf(inputs,out,0.1);
% out.T0 = an.T;
% out.FM = FM;
out.Re = Re;
out.Mach = Mach;
out.unit_system = unit_system;