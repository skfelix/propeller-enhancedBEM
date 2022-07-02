function analysis = f_adkinsPerfMorgado(inputs,geom,polares,range,beta_ref_deg,x_beta_ref)

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
%%
R = D/2;
r = x*R;
sigma = B*c./(2*pi*r);
r_hub = x_hub*R;

[rho,~,~,asound, mu] = ISA(h,unit_system);
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
else
    load('polaresNeelBright.mat','pol','foil');
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

% Initial Guesses
a = 0.1;
b = 0.1;
W = W0;

for i = length(V):-1:1
    
    J(i) = V(i)/(n*D);
    omega = n*2*pi;
    lambda = V(i)/(omega*R);
    chi = omega.*r/V(i);

    phi0 = atan(V(i)./(omega*r));   
    W0 = sqrt(V(i)^2 + (omega.*x.*R).^2);
    
    erro = 10e6;
    count = 1;
    while max(erro) > 1e-5
        
        phi = atan(V(i).*(1+a)./(omega.*r.*(1-b)));
        
        if phi(end) < 0
            phi(end) = deg2rad(1);
        end
           
        f  = 0.5*B*(1-x)./(x.*sin(phi)); % adkins
        F_tip  = abs(2/pi*acos(exp(-f)));        
        fh = 0.5*B*(x-x_hub)./(x.*sin(phi));
        F_hub = abs(2/pi*acos(exp(-fh)));        
        F  = F_tip.*F_hub;

        alpha = beta - phi;
        alphadeg = rad2deg(alpha);
        
        for j = 1:length(r)
            pol3D = postStallDuSelig(pol,R,V(i),n,r(j),c(j));
%             pol3D = postStallCorriganSchillings(pol,r(j),c(j));
            pol360 = spera360(pol3D,foil,geom,1);
%             pol360 = viternaCorrigan360(pol3D,foil);
%             plot(pol.alpha,pol.CL,pol3D.alpha,pol3D.CL,pol360.alpha,pol360.CL)
            CL(j) = interp1(pol360.alpha,pol360.CL,alphadeg(j),'linear','extrap');
            CD(j) = interp1(pol360.alpha,pol360.CD,alphadeg(j),'linear','extrap');
        end

        % CD Mach Reynolds Correction, using Re75 as reference. Hernandez &
        % Crespo (1987)
        Re = rho*W.*c/mu;
        Mach = W/asound;   
        CD = CD.*(Re75./Re).^0.2;
        CD = CD.*sqrt(1-meanMach.^2)./sqrt(1-Mach.^2);
        
        Cy = CL.*cos(phi) - CD.*sin(phi);
        Cx = CL.*sin(phi) + CD.*cos(phi);
        
        a_new = abs(1./(4*F.*sin(phi).^2./(sigma.*Cy)-1));
        b_new = 1./(4*F.*cos(phi).*sin(phi)./(sigma.*Cx)+1);
        
        
        % Radial interference velocity correction Morgado (2015)
        Wa = V(i)*(1+a_new);
        m_dot = trapz(r,2*pi*rho*Wa.*r);
        Wa_bar = m_dot/(pi*rho*R^2);
        Vt = omega*r.*b_new;
        
        Q = trapz(r,4*pi*rho*Wa_bar.*Vt.*r);
        Vt75 = 2/3*Q/(pi*rho*Wa_bar*R*(R^2-r_hub^2));
        Vtnew = 0.75*R*Vt75./r;
        b_corrected = Vtnew./(omega*r);
        b_new = b_corrected;   
  
        erro = max([abs(a-a_new); abs(b-b_new)]);
        k = 0.3;
        a = (1-k)*a + k*a_new;
        b = (1-k)*b + k*b_new;
        
        if ~isreal(a) || ~isreal(b) 
        	a = abs(a);
            b = abs(b);
        end
        
        W = sqrt((V(i).*(1+a)).^2 + (omega*r.*(1-b)).^2);
        count = count + 1;
        if count > 100
            break;
        end
    end
    
    dT = 0.5*rho.*W.^2*B.*c.*Cy;
    T(i) = trapz(r,dT)
    dQ = 0.5*rho.*W.^2*B.*c.*Cx.*r;
    Q = trapz(r,dQ);
    P = omega*Q; % 2*pi*n*Q  
    Ct(i) = T(i)./(rho*n^2*D^4);
    Cq(i) = Q./(rho*n^2*D^5);
    Cp(i) = P./(rho*n^3*D^5);
    eta(i) = J(i)*Ct(i)/Cp(i);   
end

analysis.J = J;
analysis.V = V;
analysis.Ct = Ct;
analysis.Cp = Cp;
analysis.eta = eta;
analysis.T = T;
analysis.unit_system = unit_system;

