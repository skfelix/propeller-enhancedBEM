function pol3D = postStallDuSelig(pol,R,V,n,r,c)

pol3D = pol;
omega = 2*pi*n;
x = r/R;

alpha = pol.alpha;
alpha_rad = deg2rad(alpha);
CL2D = pol.CL;
CD2D = pol.CD;
CD0 = interp1(alpha,CD2D,0);  % Cd for alpha = 0
alpha0_rad = interp1(CL2D,alpha_rad,0);
linReg = polyfit(alpha(alpha>-5 & alpha<5),CL2D(alpha>-5 & alpha<5),1);         
CLa = linReg(1); % dCL/dalpha 
CLa_rad = CLa*180/pi;

a = 1;
b = 1;
d = 1;
W = sqrt(V.^2+(omega*r)^2);
Gamma = omega*R/W;

fL = 1/(2*pi)*((1.6/0.1267*(c/r)*(a - (c/r)^(d/(Gamma*x)))/(b+(c/r)^(d/(Gamma*x))))-1);
fD = 1/(2*pi)*((1.6/0.1267*(c/r)*(a - (c/r)^(d/(2*Gamma*x)))/(b+(c/r)^(d/(2*Gamma*x))))-1);

CLp = CLa_rad*(alpha_rad - alpha0_rad); % Inviscid Thin-airfoil theory - Branlard (2011) Eq 1.25 

deltaCL = fL*(CLp - CL2D);
deltaCD = - fD*(CD2D - CD0);
CL3D = CL2D + deltaCL;
CD3D = CD2D + deltaCD;

if x >= 0.8
    CL3D = CL2D - (omega*r/W).^2*(CLp - CL2D).*CL2D./CLp;
end

pol3D.CL = CL3D;
pol3D.CD = CD3D;

% [pols] = spera360(pol,foil)
% [pol3Ds] = spera360(pol3D,foil)
% 
% figure
% hold on
% % plot(alpha,CL2D,alpha,CL3D,'--')
% plot(pols.alpha,pols.CL,pol3Ds.alpha,pol3Ds.CL,'--')
% legend('2D','3D','2D Extended','3D Extended')