function polar = spera360(polar,foil,geom,ARinf)
% [pol,foil] = xfoil('clark-y.dat',-10:0.5:20,2.3e6,0,'oper iter 150');
pol = polar;
% % function polar = spera360(polar,foil,geom,varargin)
t_c = foil.max_thickness;
% t_c = 0.12;

if ~ARinf
    c = geom.c;
    x = geom.x;
    r = x*geom.D;
    R = geom.D/2;
    
    % Blade Aspect Ratio
    Rt = R;
    Rh = geom.x_hub*R;
    Rm = sqrt(0.5*(Rt^2 + Rh^2));
    cm = interp1(x,c,Rm/R);
    
    AR = 2*(Rt-Rh)/cm;
else
    AR = inf;
end

% 2D parameters
A_prime = polar.alpha;                            % alpha
CL1_prime = polar.CL;                             % CL
A0 = interp1(CL1_prime,A_prime,0,'linear','extrap');           % alpha for CL1=0
CL1max_prime = max(CL1_prime);                              % CLmax
ACL1_prime = interp1(CL1_prime,A_prime,CL1max_prime);           % alpha for CL1
linReg = polyfit(A_prime(A_prime>-5 & A_prime<5),CL1_prime(A_prime>-5 & A_prime<5),1);         
S1_prime = linReg(1); % dCL/dalpha 
CD1_prime = polar.CD;
CDmin = min(CD1_prime);
ACDmin = interp1(CD1_prime,A_prime,CDmin);
ACD1_prime = A_prime(end);                    % Alpha for max CD
CD1max_prime = interp1(A_prime,CD1_prime,ACD1_prime); % CD at CL1


% 3D parameters - Pre-stall
E_tau = 0.316*AR^(-0.90);
E_sigma = 0.280*AR^(-0.90);

A = A_prime + pi*CL1_prime*E_tau;
CD1 = CD1_prime + CL1_prime.^2*E_sigma;
ACL1 = ACL1_prime + 18.2*CL1max_prime*AR^(-0.90);
S1 = S1_prime/(1 + 18.2*S1_prime*AR^(-0.90));
ACD1 = ACD1_prime + 18.2*CL1max_prime*AR^(-0.90);
CD1max = CD1max_prime + 0.280*CL1max_prime^2*AR^(-0.90);
CL1max = CL1max_prime*(0.67 + 0.33*exp(-(4/AR)^2));
%%
% 3D parameters - Post-stall
F1 = 1.190*(1-t_c^2);
F2 = 0.65 + 0.35*exp(-(9/AR)^2.3);
G1 = 2.300*exp(-(0.65*t_c)^0.90);
G2 = 0.52 + 0.48*exp(-(6.5/AR)^1.1);

CL2max = F1*F2;
CD2max = G1*G2;


% CL2max = 1.17*(0.67+0.33*exp(-20/AR))
% CD2max = 1.98*(0.60+0.40*exp(-20/AR))

% Curves

% Pre-stall curves
RCL1 = S1*(ACL1-A0) - CL1max;
N1 = 1+CL1max/RCL1;
RCL2 = 1.6132-CL2max;
N2 = 1 + CL2max/RCL2;

alpha = A_prime(1):0.25:120;
CL = zeros(1,length(alpha));
CD = zeros(1,length(alpha));
%
for i=1:length(alpha)
    % Lift
    if alpha(i) <= A0
        CL(i) = S1*(alpha(i)-A0)+RCL1*((A0-alpha(i))/(ACL1-A0))^N1;
    elseif alpha(i) > A0 && alpha(i) < ACD1 && AR~=inf
        CL(i) = S1*(alpha(i)-A0)-RCL1*((alpha(i)-A0)/(ACL1-A0))^N1;
    elseif alpha(i) > A0 && alpha(i) < ACD1 
        CL(i) = interp1(A_prime,CL1_prime,alpha(i));
    elseif alpha(i) >= ACD1  && alpha(i) <= 92.0
        CL(i) = -0.032*(alpha(i)-92) - RCL2*((92-alpha(i))/51)^N2;
    elseif alpha(i) > 92.0
        CL(i) = -0.032*(alpha(i)-92) + RCL2*((alpha(i)-92)/51)^N2;
    end    
    
    % Drag
%     if alpha(i) <= A_prime(1) %(2*ACDmin-ACD1)
%         CD(i) = (CD1max + (CD2max-CD1max)*sind((ACD1 - alpha(i))/(90-ACD1)*90));
%     elseif alpha(i) > (2*ACDmin-ACD1) && alpha(i) <= ACD1
    if alpha(i) >= A_prime(1) && alpha(i) <= ACD1
%         CD(i) = interp1(A,CD1,alpha(i));
            CD(i) = CDmin + (CD1max-CDmin)*((alpha(i)-ACDmin)/(ACD1-ACDmin))^4;
    elseif alpha(i) > ACD1 
        CD(i) = CD1max + (CD2max-CD1max)*sind((alpha(i)-ACD1)/(90-ACD1)*90);
    end
end


polar.alpha = alpha;
polar.CL = CL;
polar.CD = CD;

% figure(1)
% subplot(121)
% hold on
% plot(pol.alpha,pol.CL)
% plot(alpha,CL)
% 
% subplot(122)
% hold on
% plot(pol.alpha,pol.CD)
% plot(alpha,CD,'--')
% 
% abd = 5;







