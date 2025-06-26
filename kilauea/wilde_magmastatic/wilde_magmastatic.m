function bg = wilde_magmastatic(params)

% MAGMASTATIC computes magmastatic background state
    % created by Keel Wilde, May 30, 2025

% n_gas is perscribed kinematcally, H2Ofrac is specified

% This background state uses only melt density and melt viscosity from alphaMELTS
    % alphamelts lookup table accounts for solid composition but otherwise
    % arbitrary values for n_gas and XH2O (not self-consistent)

% --- gas density values based on selection ---
% Gas Equation of State (Duan & Zhang, 2006)
    % EOS of the H2O, CO2, and H2O–CO2 systems up to 10 GPa and 2573.15 K
        % interpolation range for current lookup table
        % P: 0.1 - 200 MPa
        % T: 1000 - 1500 deg C
        % XH2O: 0-1
% Ideal gas law for mixed CO2-H2O
% alphaMELTs lookup table for gas density 
    % based on 0.01 n_tot
    % NOTE: current lookup table accounts for H2Ofrac = 0.3 ONLY

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

    % NOTE: 
        % this takes input: zvec, laketop = 0, cond bot = L
        % this outputs: zvec, laketop = L, cond bot = 0

%% 

% load lookup table
params.melt_data = load('NoExsolve_H2Ofrac_1.0_1155.0C_ntot_0.0001_meltdensity_meltviscosity_only.mat');
params.gas_data = load('NoExsolve_H2Ofrac_0.3_1155.0C_ntot_0.01.mat'); % only used for alphamelts gas eos

% convert Temp celcius to temp kelvin
params.T_K = params.T_C + 273.15;

% define depth vector (full column)
zvec_col = params.z';
% we will make all depth profiles 1 x N for now and reorient them at the end
params.zvec_col = zvec_col;

%% BUILD n_gas VECTOR
% we are not using an exsolution model here, prescribe n_gas kinematically
    % right now, z=0 is lake top, z=L is conduit bottom
    % this will be flipped at the end to be consistent with Chao's model

if params.interface_split

    % in this case we will linearly interpolate between points
    % no smoothing

    % initialize gas vector
    ngas_vec = zeros(size(params.zvec_col));

    % linearly interpolate lake segment: 0 ≤ z ≤ Hlake
    idx_lake = (params.zvec_col <= params.Hlake);
    ngas_vec(idx_lake) = params.ngas_laketop ...
        + (params.ngas_lakebot - params.ngas_laketop)/params.Hlake ...
          .* params.zvec_col(idx_lake);
    
    % linearly interpolate conduit segment: Hlake < z ≤ Lcol
    idx_cond = (params.zvec_col > params.Hlake) & (params.zvec_col <= params.Lcol);
    ngas_vec(idx_cond) = params.ngas_condtop ...
        + (params.ngas_condbot - params.ngas_condtop)/(params.Lcol - params.Hlake) ...
          .* (params.zvec_col(idx_cond) - params.Hlake);
    
    params.ngas_vec = ngas_vec;

else
    
    % in this case we smooth the profile
    
    % specify anchors for transitions
    z1 = params.Hlake;   % lakebot → condtop
    z2 = params.Lcol;    % condtop → condbot
    
    % adjust coefficients for smoothing widths
    alpha1 = 0.2;  % fraction of lake depth
    alpha2 = 0.8;  % fraction of conduit length
    
    % smoothing widths (meters)
    w1 = alpha1 * params.Hlake;
    w2 = alpha2 * (params.Lcol - params.Hlake);
    
    % create smooth transition functions
    t1 = 0.1 * (1 + tanh((zvec_col - z1)/w1));  % ranges from 0 → 1
    t2 = 0.1 * (1 + tanh((zvec_col - z2)/w2));  % ranges from 0 → 1
    
    % build smooth profile
    ngas_vec = ...
        params.ngas_laketop + ...
        (params.ngas_lakebot - params.ngas_laketop) * t1 + ...
        (params.ngas_condtop - params.ngas_lakebot) * t1 + ...
        (params.ngas_condbot - params.ngas_condtop) * t2;
    
    
    params.ngas_vec = ngas_vec;
end


%% solve RHS of EOS
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol = ode45(@(t,y) RHS(t,y,params),[0 params.Lcol],[params.p_atm],options);      % solve ODE

%evaluate solution and derivatives at high order
[pvec,dPdz]=deval(zvec_col,sol); 

% re-orient pvec
% pvec = pvec';

% % initialize depth vectors
% rhovec_gas = NaN(length(pvec),1);
% rhovec_melt = NaN(length(pvec),1);
% rhovec_bulk = NaN(length(pvec),1);
% 

% extract values from lookup table

pressure_Pa_table = params.melt_data.Pressure_Pa;
melt_density_table = params.melt_data.rho_melt_kgm3; 
melt_viscosity_table = params.melt_data.melt_viscosity_Pas;
gas_density_table = params.gas_data.rho_gas_kgm3; % only used for alphamelts gas eos


% we need the original P values for interpolation
P_Pa_range = pressure_Pa_table(1:end); % Pa

% extract table values
melt_densities_lookup = melt_density_table(1:end); % kg/m3
melt_viscosities_lookup = melt_viscosity_table(1:end); % Pa s
gas_densities_lookup = gas_density_table(1:end); % kg/3 % only used for alphamelts gas eos

% Initialize output vector
rhovec_bulk = nan(size(pvec));
rhovec_melt = nan(size(pvec));
rhovec_gas = nan(size(pvec));
muvec_melt = nan(size(pvec));

% loop over the pressure values
for i = 1:length(pvec)
    
    % 1D Interpolation to find depth vectors
    % melt density 
    rhovec_melt(i) = interp1(P_Pa_range, melt_densities_lookup, pvec(i), 'linear');
    muvec_melt(i) = interp1(P_Pa_range, melt_viscosities_lookup, pvec(i), 'linear');
    
    % EOS for gas density
    switch params.gasEOS

        % ideal gas law
        case 'Idealmixed'
            
            rho_H2Ogas = pvec(i) / (params.R_H2O * params.T_K);
            rho_CO2gas = pvec(i) / (params.R_CO2 * params.T_K);
            rhovec_gas(i) = ((params.H2Ofrac/rho_H2Ogas)+((1-params.H2Ofrac)/rho_CO2gas))^(-1);

        case 'D&Zmixed'

            rhovec_gas(i) = DZ2006EOS(pvec(i),params.T_K,params.H2Ofrac);

        case 'alphamixed'

            rhovec_gas(i) = interp1(P_Pa_range, gas_densities_lookup, pvec(i), 'linear');
    end


    % calculate bulk density depth vector
    if ngas_vec(i) > 0
        rhovec_bulk(i) = ((ngas_vec(i) / rhovec_gas(i)) + ((1 - ngas_vec(i))/rhovec_melt(i)))^(-1); 
    else 
        rhovec_bulk(i) = rhovec_melt(i);
    end


end


% apparent viscosity vector (these are all shear viscosities)
porosityvec = (rhovec_melt - rhovec_bulk) ./ (rhovec_melt - rhovec_gas); % eqn 40, CK22
muvec = muvec_melt ./ (1 - porosityvec); % eqn 40, CK22

% calculate mixture soundspeed
cvec = sqrt((gradient(rhovec_bulk)./gradient(pvec)).^(-1));


% alternative calculation for bulk soundspeed (doesn't seem correct)
% K_vec = -(1./v_vec).*dvdp_vec; % mixture compressibility
% % bulk soundspeed
% cvec = (rhovec_bulk.^(-1).*(K_vec.*((-1)+ngas_vec).*pvec+(-1).*ngas_vec.*(K_vec+pvec).*params.R.*rhovec_bulk.* ...
%     params.T).^2.*((-1).*K_vec.*((-1)+ngas_vec).*pvec.^2+ngas_vec.*(K_vec+pvec).^2.*params.R.*rhovec_bulk.* ...
%     params.T).^(-1)).^(1/2);  


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
% bg.K_vec = flip(K_vec);

%% plot background state %%

% % % plot colors % % %
red = "#BF4539";
orange = "#D9863D";
green = "#BCBF65";
blue = "#51A6A6";
% % % % % % % % % % % %

figure(1)

subplot(1,5,1)
plot(bg.pvec, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('pressure Pa', "FontSize", 12);
grid

subplot(1,5,2)
plot(bg.rhovec, bg.zvec_col, 'color', red, 'LineWidth', 1.5)
hold on
plot(bg.rhovec_melt, bg.zvec_col, 'color', red,'LineStyle', '--', 'LineWidth', 1.5)
plot(bg.rhovec_gas, bg.zvec_col, 'color', red,'LineStyle', ':', 'LineWidth', 1.5)
ylabel('height in column m', "FontSize", 12);
xlabel('bulk density kg/m3', "FontSize", 12);
legend('bulk','melt', 'gas')
grid

subplot(1,5,3)
plot(bg.rhovec .* bg.cvec, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('Impedance', "FontSize", 12);
grid

subplot(1,5,4)
plot(bg.ngas_vec, bg.zvec_col, 'color', green, 'LineWidth', 1.5)
xlabel('n_{gas}', "FontSize", 12);
grid

subplot(1,5,5)
plot(bg.cvec, bg.zvec_col, 'color', blue, 'LineWidth', 1.5)
xlabel('soundspeed m/s', "FontSize", 12);
grid





end

%% RHS ODE, solve for pressure as a function of z %%

function RHS = RHS(t,y,params)

% unpack data struture
pressure_Pa_table = params.melt_data.Pressure_Pa;
melt_density_table = params.melt_data.rho_melt_kgm3; 
gas_density_table = params.gas_data.rho_gas_kgm3;

% extract table values
P_Pa_range = pressure_Pa_table(1:end);
melt_densities_lookup = melt_density_table(1:end);
gas_densities_lookup = gas_density_table(1:end); % kg/3


P=y;

% 1D Interpolation to find melt densities at P
% initialize vectors
rho_gas = nan(size(P));
rho_melt = nan(size(P));

for i = 1:length(P)
    rho_melt(i) = interp1(P_Pa_range, melt_densities_lookup, P(i), 'linear', 'extrap');
    
        % EOS for gas density
    switch params.gasEOS

        % ideal gas law
        case 'Idealmixed'
            
            rho_H2Ogas = P(i) / (params.R_H2O * params.T_K);
            rho_CO2gas = P(i) / (params.R_CO2 * params.T_K);
            rho_gas(i) = ((params.H2Ofrac/rho_H2Ogas)+((1-params.H2Ofrac)/rho_CO2gas))^(-1);

        case 'D&Zmixed'

            rho_gas(i) = DZ2006EOS(P(i),params.T_K,params.H2Ofrac);

        case 'alphamixed'
            
            rho_gas(i) = interp1(P_Pa_range, gas_densities_lookup, P(i), 'linear');
    end
    
   

end


% interpolate n_gas from the given vertical profile
n_gas = interp1(params.zvec_col, params.ngas_vec, t, 'linear', 'extrap');

% calculate bulk density
if n_gas > 0
    rho = ((n_gas ./ rho_gas) + ((1.-n_gas)./rho_melt)).^(-1); 
    else 
    rho = rho_melt;

end
% handle any potential NaNs
if any(isnan(rho))
    warning('some values of bulk density are NaN. check interpolation ranges.')
end


RHS = (rho .* params.g );


end
