
%% define the material properties and discretization parameters.
Mc.R  = .2;  % conduit radius (m)
Mc.L   = 20; % conduit length (m)
Mc.nz = 2^5;  % number of grid points in z direction

Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); %cross-sectional area.
Mc.g = 9.8; % gravitational acceleration.


Mc.mu = 1e-3*ones(Mc.nz+1, 1);

z = Mc.L/Mc.nz*[0:Mc.nz]';

Mc.z = z;

% density profile, I use exponential profile as an example but you can
% specify other profiles. You just need to make sure that this profile is
% thermodynamically stable, which means the parameter M must be
% non-negative.

rho0     = 600; %top density kg/m3
rho1     = 900; %bottom density
Rhoalpha    = (log(rho1)-log(rho0))/Mc.L;

Mc.rho = rho0*exp(Rhoalpha*(Mc.L-z)); % density

Mc.c    = 1000*ones(Mc.nz+1,1); % wave speed (constant for now) m/s

Mc.K    = Mc.rho.*Mc.c.^2; % bulk modulus
Mc.Mg = Rhoalpha - Mc.rho*Mc.g./Mc.K; % the parameter M I defined in Part I Liang et al.

%density difference and average over L
Mc.rhobar = 1/Mc.L * trapz(z(2)-z(1),Mc.rho);
Mc.deltarho = rho1-rho0;

Mc.epsilon = 1;%.01; % the area ratio between the conduit and the lava lake.

%if there is a quasistatic reservoir at the bottom, define its
%properties
%assume sphere for now
A_c = Mc.S(1);%cross sectional area of conduit
R_c = 5; %radius of spherical chamber (m)

Mc.V_c = 4/3*pi*R_c^3; %volume of chamber 

K_w = 3e9;% bulk modulus of wall rock (Pa)

nu_w = 0.25;% Poisson's ratio of wall rock

G_w = 3*K_w*(1-2*nu_w)/(2+2*nu_w);% shear modulus of rock wall rock  

K_c = 4*G_w/3; %sphere elastic stiffness 

%parameters associated with coupling - using mixture bulk modulus here
%which assumes base of the conduit is similar to chamber fluid properties

Mc.alpha = A_c/Mc.V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).

Mc.Ct = Mc.V_c *(1/Mc.K(1)+1/K_c); %elastic storativity of sphere

%% set up structure for input into the conduit-reservoir mode model
out.M = Mc;
out.z = Mc.z;

%call conduit reservoir mode reduced order model
[CRout] = CR_rom_crozierkarlstrom(out);




