function [Mc] = setparameters()

%% define the material properties and discretization parameters.
Mc.R  = 10;  % conduit radius, m
Mc.L   = 700; % conduit length, m
Mc.nz = 2^8;  % number of grid points in z direction
Mc.nr = 2^4;  % number of grid points in r direction
Mc.order = 8; % order of accuracy in z direction
Mc.order_r = 8; % order of accuracy in r direction
Mc.Rres = 1000; %reservoir radius (if spherical bc), m


Mc.T = 1155;             % magma temp (deg C)
                         % right now this is generalized for the whole system
                         % ONLY TEMP FOR LOOKUP TABLES: 1155C

Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); %cross-sectional area.
Mc.g = 9.8; % gravitational acceleration.
% I don't implement the exsolution in this model, so set this to false.
Mc.with_exsolution=false;

% This interface split flag handles the properties jump at the exsolution
% depth.  You need to specify the jump index for the code to handle the
% interface condition.
Mc.interface_split=false;


z = Mc.L/Mc.nz*[0:Mc.nz]';

Mc.z = z;


% specify background state 
% options:'parameterized', 'wilde_magmastatic'
bgstate = 'wilde_magmastatic';

Mc.epsilon = 1;%.01; % the area ratio between the conduit and the lava lake.
Mc.pT.A = 2e5; % pressure perturbation amplitude
Mc.pT.T = 1.5; % pressure perturbation duration
Mc.pT.t  = 5; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.


Mc.BCtype = 'quasistatic';%'pressure'; %set the basal boundary conditions, "pressure" is p=0

if strcmp(Mc.BCtype,'quasistatic')
    %if there is a quasistatic reservoir at the bottom, define its
    %properties
    %assume sphere for now
    A_c = Mc.S(1);%cross sectional area of conduit
    R_c = Mc.Rres; %radius of spherical chamber
    Mc.V_c = 4/3*pi*R_c^3; %volume of chamber 

    K_w = 3e9;% bulk modulus of wall rock
    nu_w = 0.25;% Poisson's ratio of wall rock
    %M.lambda_w = 3*M.K_w*M.nu_w/(1+M.nu_w); % lame constant of wall rock
    G_w = 3*K_w*(1-2*nu_w)/(2+2*nu_w);% shear modulus of rock wall rock  

    K_c = 4*G_w/3; %sphere elastic stiffness 
    
end

%%

switch bgstate


    case 'wilde_magmastatic'
        
        % prescribe total exsolved gas (H2O + CO2) kinematically
        %   specify n at top of lake, top of conduit, bottom of conduit
        params.ngas_laketop = 0.004;  % mass fraction (0.01 mass frac = 1% wt%)        
        params.ngas_condtop = 0.003;
        params.ngas_condbot = 0.00;
        
        params.Lcol = Mc.L;           % column height, m
        params.Hlake = 300;           % lake depth, m
        params.Rres = Mc.Rres;        % reservoir radius, m
        params.dz = Mc.nz;            % spatial steps, m (KW note: not sure if this variable is being used in magmastatic)
        
        params.T = Mc.T;             % magma temp (deg C)
                                     % right now this is generalized for the whole system
        
        params.H2Ofrac = 0.3;       % fraction of total volatile contents (H2O+CO2) that is H2O
        Mc.H2Ofrac = params.H2Ofrac;
        
        % boundary conditions and fixed variables
        
        params.p_atm = 1e5;                % atmospheric pressure (Pa)
        params.g = Mc.g;                   % acceleration due to gravity (m/s2)
        params.sigma = 4e-6;               % empirical constant for Henry's Law (from Kalrstrom & Dunham (2016))                                  
        params.rho0 = 3000;                                                
        params.R = 462;                    % ideal gas constant (J/kg^(-1) K)
        
        params.z = Mc.z;                   % KW note: background profiles are flipped to be compatible with z=0 bottom of conduit
        
        % solve background profile
        bg = magmastatic(params);


        Mc.rho = bg.rhovec; % density
        Mc.c    = bg.cvec; % soundspeed (m/s)
        Mc.K    = Mc.rho.*Mc.c.^2; % melt bulk mod (Pa)


        %density difference and average over L
        Mc.rhobar = 1/Mc.L * trapz(z(2)-z(1),Mc.rho);

        Mc.bg_pvec = bg.pvec; % background pressure profile in column

        Mc.drhodz = gradient(Mc.rho) ./ gradient(Mc.z);
        Mc.Mg = -(((1 ./ Mc.rho) .* (Mc.drhodz)) + ((Mc.rho .* Mc.g)./ Mc.K)); % the parameter M, defined in Chao Part I paper Eqn 15
      
        Mc.mu = bg.muvec; % shear viscosity profile, Pa
        Mc.n_gas = bg.ngas_vec; % kinematically prescribed profile of exsolved gas
        
        Mc.alpha = A_c/Mc.V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).
    
        Mc.Ct = Mc.V_c *(1/Mc.K(1)+1/K_c); %elastic storativity of sphere


    case 'parameterized'

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

        %shear viscosity profile
        Mc.mu = 100*ones(Mc.nz+1, 1);

        Mc.alpha = A_c/Mc.V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).
        
        Mc.Ct = Mc.V_c *(1/Mc.K(1)+1/K_c); %elastic storativity of sphere
        

end




end


