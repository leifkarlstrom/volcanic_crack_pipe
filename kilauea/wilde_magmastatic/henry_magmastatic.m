function bg = henry_magmastatic(params)

% MAGMASTATIC computes magmastatic background state
    % created by Keel Wilde, June 2, 2025

% n_tot and H2Ofrac are specified and gas exsolution is defined by Henry's law

% This background state uses only melt density and melt viscosity from alphaMELTS
    % alphamelts lookup table accounts for solid composition but otherwise
    % arbitrary values for n_gas and XH2O (not self-consistent)

% gas density values based on Gas Equation of State (Duan & Zhang, 2006)
    % EOS of the H2O, CO2, and H2O–CO2 systems up to 10 GPa and 2573.15 K
        % interpolation range for current lookup table
        % P: 0.1 - 200 MPa
        % T: 1000 - 1500 deg C
        % XH2O: 0-1
% OR Ideal gas law for mixed CO2-H2O


%   variables
    

    % Lcol: column height (lake + conduit)              [m]
    % Rres: reservoir radius                            [m]
    % dz: spatial steps in z direction                  [-]
    % p_atm: atmospheric pressure specified for z=0     [Pa]
    % g: acceleration due to gravity                    [m/s2]
    % zvec: depth vector; z=0 is top of the lava lake   [m]
    % pvec: pressure vector                             [Pa]
    % rhovec_bulk: bulk density vector                  [kg/m3]
    % rhovec_gas: gas density vector                    [kg/m3]
    % rhovec_melt: melt density vector                  [kg/m3]
    % sigma: empirical constant for Henry's Law         [Pa -1/2]
    % nvec_gas: volatile gas mass fraction vector       [wt%]

%% 

% load lookup table
params.melt_data = load('NoExsolve_H2Ofrac_1.0_1155.0C_ntot_0.0001_meltdensity_meltviscosity_only.mat');


% convert Temp celcius to temp kelvin
params.T_K = params.T_C + 273.15;

% define depth vector (full column)
zvec_col = params.z';
% we will make all depth profiles 1 x N for now and reorient them at the end
params.zvec_col = zvec_col;



%% solve RHS of EOS
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(@(t,y) RHS(t,y,params),[0 params.Lcol],[params.p_atm],options);      % solve ODE

%evaluate solution and derivatives at high order
[pvec,dPdz]=deval(zvec_col,sol); 


% extract values from lookup table

pressure_Pa_table = params.melt_data.Pressure_Pa;
melt_density_table = params.melt_data.rho_melt_kgm3; 
melt_viscosity_table = params.melt_data.melt_viscosity_Pas;


% we need the original P values for interpolation
P_Pa_range = pressure_Pa_table(1:end); % Pa

% extract table values
melt_densities_lookup = melt_density_table(1:end); % kg/m3
melt_viscosities_lookup = melt_viscosity_table(1:end); % Pa s


% interpolate lookup table to find rho for a given P
% Initialize output vector
rhovec_bulk = nan(size(pvec));
rhovec_melt = nan(size(pvec));
rhovec_gas = nan(size(pvec));
muvec_melt = nan(size(pvec));
ngas_vec = nan(size(pvec));

%% build profiles
    % n_tot: Henryʻs law for exsolution
    % rhovec_melt, muvec_melt: alphaMELTS lookup for melt EOS and melt viscosity
    % rhovec_gas: D&Z2006 lookup table for gas EOS OR Ideal gas law
        

for i = 1:length(pvec)

    % 1D Interpolation to find depth vectors
    rhovec_melt(i) = interp1(P_Pa_range, melt_densities_lookup, pvec(i), 'linear');
    muvec_melt(i) = interp1(P_Pa_range, melt_viscosities_lookup, pvec(i), 'linear');


    if pvec(i) < (params.n_tot / params.sigma)^2
        ngas_vec(i) = (params.n_tot-params.sigma .*pvec(i).^(1/2));

        % EOS for gas density
        switch params.gasEOS
    
            % ideal gas law
            case 'Idealmixed'
                
                rho_H2Ogas = pvec(i) / (params.R_H2O * params.T_K);
                rho_CO2gas = pvec(i) / (params.R_CO2 * params.T_K);
                rhovec_gas(i) = ((params.H2Ofrac/rho_H2Ogas)+((1-params.H2Ofrac)/rho_CO2gas))^(-1);
    
            case 'D&Zmixed'
    
                rhovec_gas(i) = DZ2006EOS(pvec(i),params.T_K,params.H2Ofrac);
        end
        
        % bulk density
        % display("P: " + pvec(i));
        % display("ngas_vec: " + ngas_vec(i));
        % display("rho_H2Ogas: " + rho_H2Ogas);
        % display("rho_CO2gas: " + rho_CO2gas);
        % display("rhovec_gas: " + rhovec_gas(i));
        % display("rhovec_melt: " + rhovec_melt(i));
        rhovec_bulk(i) = ((ngas_vec(i) / rhovec_gas(i)) + ((1 - ngas_vec(i))/rhovec_melt(i)))^(-1); 
        % display("rhovec_bulk: " + rhovec_bulk(i));
        % keyboard
    else 
        % no gas exsolved
        ngas_vec(i) = 0;
        disp("No gas exsolved at " + zvec_col(i))
    
        % so no gas density
        rhovec_gas(i) = 0;
        % so bulk density == melt density here
        rhovec_bulk(i) = rhovec_melt(i);
    end
end



% apparent viscosity vector (these are all shear viscosities)
porosityvec = (rhovec_melt - rhovec_bulk) ./ (rhovec_melt - rhovec_gas); % eqn 40, CK22
muvec = muvec_melt ./ (1 - porosityvec); % eqn 40, CK22

% calculate mixture soundspeed
cvec = sqrt((gradient(rhovec_bulk)./gradient(pvec)).^(-1));

% save background state
% in Chao's model, z=0 is bottom of column, so we flip all of these
% profiles here
bg.zvec_col = zvec_col';
bg.pvec = flip(pvec)';
bg.rhovec = flip(rhovec_bulk)';
bg.muvec = flip(muvec)';
bg.cvec = flip(cvec)';
bg.rhovec_melt = flip(rhovec_melt)';
bg.rhovec_gas = flip(rhovec_gas)';
bg.ngas_vec = flip(ngas_vec)';

%% plot background state %%

% % % plot colors % % %
red = "#BF4539";
orange = "#D9863D";
green = "#BCBF65";
blue = "#51A6A6";
% % % % % % % % % % % %

figure(1)
subplot(1,5,1)
plot(bg.rhovec, bg.zvec_col, 'color', red, 'LineWidth', 1.5)
ylabel('height in column m', "FontSize", 12);
xlabel('bulk density kg/m3', "FontSize", 12);
grid
% set(gca, 'YDir','reverse')
% 
subplot(1,5,2)
plot(bg.pvec, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('pressure Pa', "FontSize", 12);
grid

subplot(1,5,3)
plot(bg.rhovec .* bg.cvec, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('Impedance', "FontSize", 12);
grid

subplot(1,5,4)
plot(bg.ngas_vec, bg.zvec_col, 'color', green, 'LineWidth', 1.5)
xlabel('n_{gas}', "FontSize", 12);
xlim([0 0.02]);
grid

subplot(1,5,5)

plot(bg.cvec, bg.zvec_col, 'color', blue, 'LineWidth', 1.5)
% xlim([0 400])

xlabel('soundspeed m/s', "FontSize", 12);
grid


figure(2)
subplot(1,2,1)
plot(bg.rhovec_melt, bg.zvec_col, 'color', red, 'LineWidth', 1.5)
ylabel('height in column m', "FontSize", 12);
xlabel('melt density kg/m3', "FontSize", 12);
ylim([0 bg.zvec_col(end)])
grid
% set(gca, 'YDir','reverse')
% 
subplot(1,2,2)
plot(bg.rhovec_gas, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('gas density kg/m3', "FontSize", 12);
ylim([0 bg.zvec_col(end)])
grid

keyboard

end

%% RHS ODE, solve for pressure as a function of z %%

function RHS = RHS(t,y,params)

P=y;

% unpack data struture for melt densities (from alphamelts)
pressure_Pa_table = params.melt_data.Pressure_Pa;
melt_density_table = params.melt_data.rho_melt_kgm3; 

% extract table values
P_Pa_range = pressure_Pa_table(1:end);
melt_densities_lookup = melt_density_table(1:end);

% Henry's Law
if P(1) < (params.n_tot / params.sigma)^2
    n_gas = (params.n_tot-params.sigma .* P.^(1/2)); 

    % EOS for gas density
    switch params.gasEOS

        % ideal gas law
        case 'Idealmixed'
            rho_H2Ogas = P ./ (params.R_H2O * params.T_K);
            rho_CO2gas = P ./ (params.R_CO2 * params.T_K);
            rho_gas = ((params.H2Ofrac/rho_H2Ogas)+((1-params.H2Ofrac)/rho_CO2gas))^(-1);

        case 'D&Zmixed'

            rho_gas = DZ2006EOS(P,params.T_K,params.H2Ofrac);
    end

else 
     % no gas exsolved
    n_gas = 0;

    % so no gas density
    rho_gas = 0;

end

% melt density
% rho_melt = params.rho0.*(1+ P./params.K_melt); % if you use this you need to update melt bulk mod
rho_melt = interp1(P_Pa_range, melt_densities_lookup, P, 'linear', 'extrap');

%bulk density
if n_gas > 0
    rho = ((n_gas ./ rho_gas) + ((1.-n_gas)./rho_melt)).^(-1); 

else 
    rho = rho_melt;

end



if any(isnan(rho))
    warning('Some requested values are outside the interpolation range.');
end
   

RHS = (rho.*params.g );


end
