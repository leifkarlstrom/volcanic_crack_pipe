% conduit internal gravity model driver
clear
close all

%% Model parameters.
% conduit parameters
%
source = '../source';
addpath(genpath(source));

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

Mc.mu = 100*ones(Mc.nz+1, 1);

z = Mc.L/Mc.nz*[0:Mc.nz]';

% density profile, I use exponential profile as an example but you can
% specify other profiles. You just need to make sure that this profile is
% thermodynamically stable, which means the parameter M must be
% non-negative.

rho0     = 1000;
rho1     = 2000;
Rhoalpha    = (log(rho1)-log(rho0))/Mc.L;

Mc.rho = rho0*exp(Rhoalpha*(Mc.L-z)); % density
Mc.c    = 1000*ones(Mc.nz+1,1); % wave speed
Mc.K    = Mc.rho.*Mc.c.^2; % bulk modulus
Mc.Mg = Rhoalpha - Mc.rho*Mc.g./Mc.K; % the parameter M I defined in Part I paper.

%density difference and average over L
Mc.rhobar = 1/Mc.L * trapz(z(2)-z(1),Mc.rho);
Mc.deltarho = rho1-rho0;


% Here I show you how to specify a jump in properties, let's assume the
% density increases from 800 to 1500 exponentially and jump to 2000 at the
% middle point and then increase to 3000 exponentially.

% More complext background properties can be used as long you specify the
% properties and the split_index.

% You also need to make sure each segment has enough grid points for
% specifed accuracy (you can increase the total number of grid points.)


% Make sure all the material properties have the same dimension
% S, mu rho, c, K, Mg.
%

if Mc.interface_split
    
    rho0      = 800;
    rho1      = 1500;
    rho2      = 2000;
    rho3      = 3000;
    
    c0         = 1000; % wavespeed in upper section
    c1         = 2000; % wavespeed in lower section
    
    mu0     = 100;% viscosity in upper section
    mu1     = 50;% viscosity in upper section
    
    z_split        = 500; % the z-coordinate of jump point.
    % find the closest point to the specified coordinate and duplicate the grid point
    [~, Mc.split_index] = min(abs(z - z_split));
    z_upper = z(Mc.split_index:end);
    z_lower = z(1:Mc.split_index);
    
    %
    L_lower = z(Mc.split_index);
    L_upper = Mc.L - L_lower;
    
    Rhoalpha_upper =  (log(rho1)-log(rho0))/L_upper*ones(size(z_upper));
    Rhoalpha_lower =  (log(rho3)-log(rho2))/L_lower*ones(size(z_lower));
    
    rho_upper = rho0*exp(Rhoalpha_upper.*(Mc.L-z_upper)); % density upper segment
    rho_lower =  rho2*exp(Rhoalpha_lower.*(L_lower-z_lower)); % density lower segment
    
    c_upper = c0*ones(size(z_upper));
    c_lower = c1*ones(size(z_lower));
    
    mu_upper = mu0*ones(size(z_upper));
    mu_lower = mu1*ones(size(z_lower));
    
    Mc.rho = [rho_lower; rho_upper];
    Mc.c    = [c_lower; c_upper];
    Mc.mu = [mu_lower; mu_upper];
    Mc.K    = Mc.rho.*Mc.c.^2; % bulk modulus
    Rhoalpha = [Rhoalpha_lower;Rhoalpha_upper];
    Mc.Mg = Rhoalpha - Mc.rho*Mc.g./Mc.K; % the parameter M I defined in Part I paper.
    z = [z_lower;z_upper];
    Mc.S = pi*Mc.R^2*ones(length(z), 1); %cross-sectional area.
end

% plot material properties

plot_properties = 0;

if plot_properties
    figure(111);
    props = {'c','rho','Mg'};
    n_props = length(props);
    for i = 1: n_props
        prop = props{i};
        subplot(1, n_props, i);
        if strcmp(prop, 'Mg')
            plot(log10(Mc.(prop)), z);
            xlabel('Log10(Mg)');
        else
            plot(Mc.(prop), z);
            xlabel(prop);
        end
        ylabel('z');
    end    
end

%%
Mc.epsilon = 1;%.01; % the area ratio between the conduit and the lava lake.
Mc.pT.A = 2e4; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t  = 5; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.

Mc.BCtype = 'quasistatic';%'pressure'; %set the basal boundary conditions, "pressure" is p=0

if strcmp(Mc.BCtype,'quasistatic')
    %if there is a quasistatic reservoir at the bottom, define its
    %properties
    %assume sphere for now
    A_c = Mc.S(1);%cross sectional area of conduit
    R_c = 750; %radius of chamber
    V_c = 4/3*pi*R_c^3; %volume of chamber 

    K_w = 3e9;% bulk modulus of wall rock
    nu_w = 0.25;% Poisson's ratio of wall rock
    %M.lambda_w = 3*M.K_w*M.nu_w/(1+M.nu_w); % lame constant of wall rock
    G_w = 3*K_w*(1-2*nu_w)/(2+2*nu_w);% shear modulus of rock wall rock  

    K_c = 4*G_w/3; %sphere elastic stiffness 
    Mc.alpha = A_c/V_c*1/((1/Mc.K(1)+1/K_c));% coupling parameters, dp_c/dt=alpha*v(0).
end

Model = conduit_internal_g(Mc);
%% time domain simulation.
CFL = 0.5;
skip = 10;
T = 80;
use_imex = true; % a flag to choose if to use IMEX or not.
plot_simu = false;

% time stepping
hmin = min([Model.geom.dz]);
cmax = max(Model.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);
time  = [skip:skip:nt]*dt;

if ~ use_imex
    A = Model.Ae + Model.Ai;
    
    % this matrix A is what you need for analyzing the eigenmode and
    % energetics.
    
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end

fun = @(u,t) A*u + Model.Fp(:,1)*Model.M.G(t);
tic

% Storage arrays
%out.p = zeros((Mc.nz+1),nt+1);
%ICs
out.p(:,1) = zeros((Mc.nz+1),1);
out.t(1) = 0; 
out.z = z;

for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
    end
    
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    
    if mod(i,skip) == 0
        iter = i/skip;
        %plot solution.
        nr  = Model.geom.nr;
        nz = Model.geom.nz;
        
        % velocity
        vz  = reshape(Model.field(Model.u, 1), [nr,nz]); % [nr, nz]
        
        % pressure
        pz  = Model.field(Model.u, 2); % [nz, 1]
        
        % displacement
        h    = Model.field(Model.u, 3); % [nz,  1]

        % lake height
        hL    = Model.field(Model.u, 4); % [nz,  1]
        
        % width averaged velocity.
        uz  = (Model.op.W1*vz)';%   [nz, 1]

        if strcmp(Mc.BCtype,'quasistatic')
            p_c = Model.field(Model.u, 5);

            out.p_c(:,iter+1) = p_c;
        end
        
        z   = Model.geom.z;
        rm  = Model.geom.rm;
        
        time_i = i*dt;
        if plot_simu
            figure(2);
            subplot(1,3,1)
            pcolor(rm, z, vz');
            shading INTERP;
            cmap;
            caxis([-1 1]);
            colorbar
            xlabel('r (m)'), ylabel('z (m)')
            
            subplot(1,3,2)
            plot(uz, z);
            xlim([-1, 1])
            ylim([0 Model.M.L])
            xlabel('u'), ylabel('z (m)')
            
            subplot(1,3,3)
            plot(pz, z);
            xlim([-2.5e5,2.5e5]);
            ylim([0 Model.M.L])
            xlabel('p'), ylabel('z (m)')
            drawnow;
        end

        out.p(:,iter+1) = pz;
        out.t = [out.t, time_i];
    end
end
