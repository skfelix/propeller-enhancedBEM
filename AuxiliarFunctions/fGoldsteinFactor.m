function [ G , vr] = fGoldsteinFactor(l_bar ,B,N)
% Computes Goldstein factor G
% l_bar: dimensionless torsional parameter h/(2 piR)
% B: number of blades
% N: number of radial positions
% AUTHOR: E. Branlard
R = 1 ; w = 1; % Using dimensionless quantities
vr = linspace (0,1,N); % Radial positions
G = fGoldsteinCirculation(l_bar*R,w,R,B,vr)*B/(2*pi*l_bar*R*w);
end % function fGoldsteinFactor
function [ GammaGoldstein ] = fGoldsteinCirculation(l_bar ,w,R,B,vr)
% Computes Goldstein circulation using superposition of helix
% AUTHOR: E. Branlard
vr = vr(2:end); % skipping point 0
N = length(vr);
vrCP = (3/(2*N):1/N:1)*R; % Control Points (CPs) in between vortices
% Calculation of matrices of influence at CPs for a helical filament
A=zeros(N,N); %
Gamma_h = 1; % Circulation is unitary
psi_h = 0; % Azimutal position is 0 (on the Lifting line)
for j=1: length(vr) % loop on helical vortex radii
for i=1: length(vrCP) % loop on Control points radii
A(i,j)=fUi_HelixNTheory(Gamma_h ,vrCP(i),vr(j),l_bar*R,psi_h ,B);
end
end
% Boundary conditions values on the vortex sheet CP
U=zeros(N,1);
U(1:(N-1))= w*1./(1+ l_bar ^2./( vrCP (1:(N-1))/R).^2);
A(end ,:) =1; U(end)=0; % Condition: sum of gamma =0
%Solving for individual trailed vorticity and cumulative one
Gamma_t=A\U; GammaGoldstein=cumsum(Gamma_t); GammaGoldstein =[0;
GammaGoldstein ];
end % function fGoldsteinCirculation
function [ uz ] = fUi_HelixNTheory(Gamma ,r,r0 ,l,psih ,B);
% Computes induction from N-helices
% AUTHOR: E. Branlard
C0z= (l^2+r0^2) ^(1/4) /(l^2+r^2) ^(1/4);
C1z= l/24*( (3*r^2-2*l^2)/(l^2+r^2) ^(3/2) +(2*l^2+9* r0^2)/(l^2+r0^2) ^(3/2) );
pexi=r/r0*(l+sqrt(l^2+r0^2))/(l+sqrt(l^2+r^2))*exp(sqrt(l^2+r^2)/l)/exp(sqrt(l^2+r0^2)/l);
mexi =1/ pexi; t=psih;
if(abs(r)<r0)
tmp = 1/( (mexi*exp(-1i*t))^B -1 );
vz = 1/(2* pi*l) + 1/(2* pi*l)*C0z*real( tmp + C1z/B*log (1+tmp) );
elseif(abs(r)>r0)
tmp = 1/( (pexi*exp(-1i*t))^B -1 );
vz = 0 + 1/(2* pi*l)*C0z*real( - tmp + C1z/B*log (1+tmp) );
else
vz=0;
end
uz = - B*Gamma*vz;
end % function fUi_HelixNTheory