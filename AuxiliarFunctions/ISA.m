% International Standard Atmosphere (ISA) in SI units
% [rho,T,p,a,mu,nu,AD]=ISA(1000, Toffset,'metric')
% Input : h = altitute (m), T_offset = temperature deviation (ºC), unit_system = 'metric' standard or 'imperial'
% Outputs : T = temperature (K), p = pressure (N/mˆ2), rho = density (kg/mˆ3), a = sound speed (m/s)

function [rho,T,p,a,mu,nu,AD]=ISA(h, varargin)

if nargin == 3
    T_offset = varargin{1};
    unit_system = varargin{2};
    if strcmp(unit_system,'imperial')
        h = h*0.3048;
    end
elseif nargin == 2
    T_offset = varargin{1};
    unit_system = 'metric';
else
    T_offset = 0;
    unit_system = 'metric';
end

g = 9.80665;
h = h/1000; % m -> km

h1 = 11;
h2 = 20;
h3 = 32;

L0 = -6.5e-3; 
L2 = 1e-3; 

g0 = 9.80665;
m0 = 28.96442;
R0 = 8314.32; 
R = R0/m0;

T0 = 288.15;
p0 = 1.01325e5;
rho0 = 1.2250;

T1 = T0+L0*h1*1e3;
p1 = p0*(T1/T0)^(-g0/(R*L0)); 
rho1 = rho0*(T1/T0)^(-(1+g0/(R*L0))); 

T2 = T1;
p2 = p1*exp(-g0/(R*T2)*(h2-h1)*1e3);
rho2 = rho1*exp(-g0/(R*T2)*(h2-h1)*1e3);

if h<=h1
	% Troposphere:
	T_ISA = T0+L0*h*1e3;
	p_ISA = p0*(T_ISA/T0)^(-g0/(R*L0));
    rho_ISA = rho0*(T_ISA/T0)^(-(1+g0/(R*L0)));
elseif h<=h2
	% Tropopause and low stratosphere:
	T_ISA = T1;
	p_ISA = p1*exp(-g0/(R*T_ISA)*(h-h1)*1e3);
	rho_ISA = rho1*exp(-g0/(R*T_ISA)*(h-h1)*1e3);
elseif h<=h3
	% Stratosphere:
	T_ISA = T2+L2*(h-h2)*1e3;
	p_ISA = p2*(T_ISA/T2)^(-g0/(R*L2));
	rho_ISA = rho2*(T_ISA/T2)^(-(1+g0/(R*L2)));
end

% Temperature offset


T = T_ISA + T_offset;
rho = rho_ISA*(T_ISA/T); % from p1 = (rho*R*T)1 = p2 = (rho*R*T)2
p = p_ISA;

if h<=h1
	% Troposphere:
    TAD = T0*(rho/rho0)^(1/(-(1+g0/(R*L0))));
    AD = (TAD-T0)/L0;
elseif h> h1
	fprintf('Density altitude can not be calculated for h > 11 0000 m \n')
    AD = NaN;
end


% Sound speed
gamma = 1.4;
a = sqrt(gamma*R*T);

% Viscosity

% Dry air - for 250 <= T <= 600 K
MA0 = -9.8601e-1;
MA1 =  9.080125e-2;
MA2 = -1.17635575e-4;
MA3 =  1.2349703e-7;
MA4 = -5.7971299e-11;
mu = (MA0 + MA1*T + MA2*T.^2 + MA3*T.^3 + MA4*T.^4)*1e-6; % [N s m^-2]
nu = mu/rho;


%% Output conversion

if nargin > 1 && strcmp(unit_system,'imperial')
    lb2kg = 0.4536;
    slug2kg = 14.5939;
    ft2m  = 0.3048;
    rho = rho*(1/slug2kg)/(1/ft2m^3); % Kg/m^3 to slug/ft^3 
    p = p*(1/(lb2kg*g))/(1/ft2m)^2; % Pa(N/m^2, kg*g/m^2) to (lbf/ft^2)
    a = a*(1/ft2m); % m/s to ft/s
    mu = mu*(1/(lb2kg*g))/(1/ft2m)^2; % Pa*s
end
