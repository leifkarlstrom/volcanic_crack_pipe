function [Mc] = setparameters(bgstate)

%% define the material properties and discretization parameters.
Mc.R  = 5;  % conduit radius
Mc.L   = 600; % conduit length
Mc.nz = 2^7;  % number of grid points in z direction
Mc.nr = 2^4;  % number of grid points in r direction
Mc.order = 6; % order of accuracy in z direction
Mc.order_r = 6; % order of accuracy in r direction

Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); %cross-sectional area.
Mc.g = 9.8; % gravitational acceleration.

% I don't implement the exsolution in this model, so set this to false.
Mc.with_exsolution=false;

% This interface split flag handles the properties jump at the exsolution
% depth.  You need to specify the jump index for the code to handle the
% interface condition.
Mc.interface_split=false;

Mc.mu = 50*ones(Mc.nz+1, 1);

z = Mc.L/Mc.nz*[0:Mc.nz]';

Mc.z = z;
% density profile, I use exponential profile as an example but you can
% specify other profiles. You just need to make sure that this profile is
% thermodynamically stable, which means the parameter M must be
% non-negative.

rho0     = 500;
rho1     = 2200;
Rhoalpha    = (log(rho1)-log(rho0))/Mc.L;

Mc.rho = rho0*exp(Rhoalpha*(Mc.L-z)); % density
Mc.c    = 1000*ones(Mc.nz+1,1); % wave speed
Mc.K    = Mc.rho.*Mc.c.^2; % bulk modulus
Mc.Mg = Rhoalpha - Mc.rho*Mc.g./Mc.K; % the parameter M I defined in Part I paper.

%density difference and average over L
Mc.rhobar = 1/Mc.L * trapz(z(2)-z(1),Mc.rho);
Mc.deltarho = rho1-rho0;

Mc.epsilon = 1;%.01; % the area ratio between the conduit and the lava lake.
Mc.pT.A = 2e4; % pressure perturbation amplitude
Mc.pT.T = 2.5; % pressure perturbation duration
Mc.pT.t  = 5; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.

Mc.BCtype = 'quasistatic';%'pressure'; %set the basal boundary conditions, "pressure" is p=0

if strcmp(Mc.BCtype,'quasistatic')
    %if there is a quasistatic reservoir at the bottom, define its
    %properties
    %assume sphere for now
    A_c = Mc.S(1);%cross sectional area of conduit
    R_c = 1000; %radius of spherical chamber
    Mc.V_c = 4/3*pi*R_c^3; %volume of chamber 

    K_w = 3e9;% bulk modulus of wall rock
    nu_w = 0.25;% Poisson's ratio of wall rock
    %M.lambda_w = 3*M.K_w*M.nu_w/(1+M.nu_w); % lame constant of wall rock
    G_w = 3*K_w*(1-2*nu_w)/(2+2*nu_w);% shear modulus of rock wall rock  

    K_c = 4*G_w/3; %sphere elastic stiffness 
    Mc.alpha = A_c/Mc.V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).
    
    Mc.Ct = Mc.V_c *(1/Mc.K(1)+1/K_c); %elastic storativity of sphere
end


