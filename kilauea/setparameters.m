function [Mc] = setparameters()

%% define the material properties and discretization parameters.

Mc.R  = 5;  % conduit radius
Mc.L   = 700; % conduit length
Mc.nz = 2^8;  % number of grid points in z direction
Mc.nr = 2^4;  % number of grid points in r direction
Mc.order = 8; % order of accuracy in z direction
Mc.order_r = 8; % order of accuracy in r direction

Mc.Rres = 1000; %reservoir radius (if spherical bc), m


Mc.T_C = 1155;             % magma temp (deg C)
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
% options:'parameterized', 'wilde_magmastatic', 'henrys law', 'alphaEOS' 
bgstate = 'wilde_magmastatic';

Mc.epsilon = 1;%.01; % the area ratio between the conduit and the lava lake.
Mc.pT.A = 2e5; % pressure perturbation amplitude %2e5
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
        % choose gas EOS
            % D&Z2006 gas density values: 'D&Zmixed'
            % mixed ideal gas law  gas density values: 'Idealmixed'
            % alphaMELTs lookup table for gas density: 'alphamixed'
                % NOTE: current lookup table accounts for H2Ofrac = 0.3 ONLY
        
        % prescribe total exsolved gas (H2O + CO2) kinematically
        %   specify n at top of lake, top of conduit, bottom of conduit
        params.ngas_laketop = 0.001;  % mass fraction (0.01 mass frac = 1% wt%) 
        params.ngas_lakebot = 0.001;
        params.ngas_condtop = 0.001;
        params.ngas_condbot = 0.0005;
        
        params.Lcol = Mc.L;           % column height, m
        params.Hlake = 300;           % lake depth, m
        Mc.Hlake = params.Hlake; 
        params.Rres = Mc.Rres;        % reservoir radius, m
        params.dz = Mc.nz;            % spatial steps, m (KW note: not sure if this variable is being used in magmastatic)
        
        params.T_C = Mc.T_C;             % magma temp (deg C)
                                     % right now this is generalized for the whole system
    
        
        params.gasEOS = 'Idealmixed';                             
        params.H2Ofrac = 0.3;       % fraction of total volatile contents (H2O+CO2) that is H2O
        Mc.H2Ofrac = params.H2Ofrac;
        
        % boundary conditions and fixed variables
        
        params.p_atm = 1e5;                % atmospheric pressure (Pa)
        params.g = Mc.g;                   % acceleration due to gravity (m/s2)
        params.sigma = 4e-6;               % empirical constant for Henry's Law (from Kalrstrom & Dunham (2016))                                  
        params.R_H2O = 462;                    % ideal gas constant water vapor (J/kg^(-1) K)
        params.R_CO2 = 189;                    % ideal gas constant CO2 vapor (J/kg^(-1) K)
        %params.rho0 = 3000;                                                
        
        params.z = Mc.z;                   % KW note: background profiles are flipped to be compatible with z=0 bottom of conduit
        
        % solve background profile
        bg = wilde_magmastatic(params);
      

        Mc.rho = bg.rhovec; % density
        Mc.c    = bg.cvec; % soundspeed (m/s)
        Mc.K    = Mc.rho.*Mc.c.^2; % melt bulk mod (Pa)
        
        %save out ngas intervals
        Mc.ngas_laketop = params.ngas_laketop; 
        Mc.ngas_condtop = params.ngas_condtop; 
        Mc.ngas_condbot = params.ngas_condbot;

        %density difference and average over L
        Mc.rhobar = 1/Mc.L * trapz(z(2)-z(1),Mc.rho);

        Mc.bg_pvec = bg.pvec; % background pressure profile in column

        Mc.drhodz = gradient(Mc.rho) ./ gradient(Mc.z);
        Mc.Mg = -(((1 ./ Mc.rho) .* (Mc.drhodz)) + ((Mc.rho .* Mc.g)./ Mc.K)); % the parameter M, defined in Chao Part I paper Eqn 15
      
        Mc.mu = bg.muvec; % shear viscosity profile, Pa
        Mc.n_gas = bg.ngas_vec; % kinematically prescribed profile of exsolved gas
        
        Mc.alpha = A_c/Mc.V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).
    
        Mc.Ct = Mc.V_c *(1/Mc.K(1)+1/K_c); %elastic storativity of sphere

    case 'henrys law'
        % uses henry's law for exsolution
        % uses alphaMELTs for melt density and melt viscosity
        
        % choose gas EOS
            % uses D&Z2006 gas density values: 'D&Zmixed'
            % uses mixed ideal gas law  gas density values: 'Idealmixed'
        params.gasEOS = 'D&Zmixed';
  
        params.n_tot = 0.005;  % mass fraction TOTAL volatiles (0.01 mass frac = 1% wt%) 
        params.H2Ofrac = 0.3;       % fraction of total volatile contents (H2O+CO2) that is H2O
        Mc.H2Ofrac = params.H2Ofrac;
        
        params.Lcol = Mc.L;           % column height, m
        % params.Hlake = 300;           % lake depth, m
        % Mc.Hlake = params.Hlake; 
        params.Rres = Mc.Rres;        % reservoir radius, m
        params.dz = Mc.nz;            % spatial steps, m (KW note: not sure if this variable is being used in magmastatic)
        
        params.T_C = Mc.T_C;             % magma temp (deg C)
                                     % right now this is generalized for the whole system
       
        
        % boundary conditions and fixed variables
        
        params.p_atm = 1e5;                % atmospheric pressure (Pa)
        params.g = Mc.g;                   % acceleration due to gravity (m/s2)
        params.sigma = 4e-6;               % empirical constant for Henry's Law (from Kalrstrom & Dunham (2016))                                                                                 
        params.R_H2O = 462;                    % ideal gas constant water vapor (J/kg^(-1) K)
        params.R_CO2 = 189;                    % ideal gas constant CO2 vapor (J/kg^(-1) K)
        
        params.z = Mc.z;                   % KW note: background profiles are flipped to be compatible with z=0 bottom of conduit
        
        % solve background profile
        bg = henry_magmastatic(params);
      

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


            if Mc.interface_split
                % CURRENTLY ALMOST DIRECTLY COPIED FROM CONDUIT_INTERNALG,
                % NEEDS TO BE MODIFIED TO MAKE CONSISTENT WITH THE NEW BG
                % STATE - LK

                rho0      = 800;
                rho1      = 1500;
                rho2      = 2000;
                rho3      = 3000;
                
                c0         = 1000; % wavespeed in upper section
                c1         = 2000; % wavespeed in lower section
                
                mu0     = 100;% viscosity in upper section
                mu1     = 50;% viscosity in lwer section

                R0 = 100; %radius in upper section
                R1 = 10; %radius in lower section
                
                z_split        = 500; % the z-coordinate of jump point.
                % find the closest point to the specified coordinate and duplicate the grid point
                [~, Mc.split_index] = min(abs(z - z_split));
                z_upper = z(Mc.split_index:end);
                z_lower = z(1:Mc.split_index);
                
                % define fluid properties and geometry for each section
                L_lower = z(Mc.split_index);
                L_upper = Mc.L - L_lower;

                S_lower = pi*R0^2*ones(size(z_lower));
                S_upper = pi*R1^2*ones(size(z_upper));
                
                Rhoalpha_upper =  (log(rho1)-log(rho0))/L_upper*ones(size(z_upper));
                Rhoalpha_lower =  (log(rho3)-log(rho2))/L_lower*ones(size(z_lower));
                
                rho_upper = rho0*exp(Rhoalpha_upper.*(Mc.L-z_upper)); % density upper segment
                rho_lower =  rho2*exp(Rhoalpha_lower.*(L_lower-z_lower)); % density lower segment
                
                c_upper = c0*ones(size(z_upper));
                c_lower = c1*ones(size(z_lower));
                
                mu_upper = mu0*ones(size(z_upper));
                mu_lower = mu1*ones(size(z_lower));
                
                %now combine to make BG grid with two sections
                Mc.rho = [rho_lower; rho_upper];
                Mc.c    = [c_lower; c_upper];
                Mc.mu = [mu_lower; mu_upper];
                Mc.K    = Mc.rho.*Mc.c.^2; % bulk modulus
                Rhoalpha = [Rhoalpha_lower;Rhoalpha_upper];
                Mc.Mg = Rhoalpha - Mc.rho*Mc.g./Mc.K; % the parameter M I defined in Part I paper.
                z = [z_lower;z_upper];
                Mc.S = [S_lower;S_upper];%pi*Mc.R^2*ones(length(z), 1); %cross-sectional area.
            end

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


