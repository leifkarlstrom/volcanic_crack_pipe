function bg = magmastatic(params)

% MAGMASTATIC computes magmastatic background state
    % created by Keel Wilde, Feb. 28, 2024
    % updated by Keel Wilde, Feb. 28, 2025

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
data = load('NoExsolve_H2Ofrac_0.3_1155.0C_ntot_0.005.mat');

% convert Temp celcius to temp kelvin
params.T = params.T + 273.15;

% define depth vector (full column)
zvec_col = params.z;
params.zvec_col = zvec_col;

%% BUILD n_gas VECTOR
% we are not using an exsolution model here, prescribe n_gas kinematically
    % right now, z=0 is lake top, z=L is condbot 
    % this will be flipped at the end to be consistent with Chao's model

% z values consistent with n_gas piecewise values
z_known = [0, params.Hlake, params.Lcol]; 
% corresponding ngas values
n_known = [params.ngas_laketop, params.ngas_condtop, params.ngas_condbot];


% interpolate to get a smooth piecewise ngas profile
% in the future, it might be better to do a linear interp and then smooth
% the transition at Hlake
ngas_vec = interp1(z_known, n_known, zvec_col, 'pchip');

params.ngas_vec = ngas_vec;


%% solve RHS of EOS
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(@(t,y) RHS(t,y,params, data),[0 params.Lcol],[params.p_atm],options);      % solve ODE

%evaluate solution and derivatives at high order
[pvec,dPdz]=deval(zvec_col,sol); 

% re-orient pvec
pvec = pvec';

% % initialize depth vectors
% rhovec_gas = NaN(length(pvec),1);
% rhovec_melt = NaN(length(pvec),1);
% rhovec_bulk = NaN(length(pvec),1);
% 

% extract values from lookup table

pressure_Pa_table = data.Pressure_Pa;
melt_density_table = data.rho_melt_kgm3; 
gas_density_table = data.rho_gas_kgm3; 
melt_viscosity_table = data.melt_viscosity_Pas;
v_table = data.bulk_volume_cm3;
dvdp_table = data.dvdp;

% we need the original P values for interpolation
P_Pa_range = pressure_Pa_table(1:end); % Pa

% extract table values
melt_densities_lookup = melt_density_table(1:end); % kg/m3
gas_densities_lookup = gas_density_table(1:end); % kg/3
melt_viscosities_lookup = melt_viscosity_table(1:end); % Pa s

% NOTE: this v and dvdp is taken from alphamelts and is not
% consistent with this n_gas profile
v_lookup = v_table(1:end); % cm3
dvdp_lookup = dvdp_table(1:end);


% interpolate lookup table to find rho for a given P
% Initialize output vector
rhovec_bulk = nan(size(pvec));
rhovec_melt = nan(size(pvec));
rhovec_gas = nan(size(pvec));
muvec_melt = nan(size(pvec));
v_vec = nan(size(pvec));
dvdp_vec = nan(size(pvec));

% loop over the pressure values
for i = 1:length(pvec)
    % 1D Interpolation to find depth vectors
    rhovec_melt(i) = interp1(P_Pa_range, melt_densities_lookup, pvec(i), 'linear');
    rhovec_gas(i) = interp1(P_Pa_range, gas_densities_lookup, pvec(i), 'linear');
    muvec_melt(i) = interp1(P_Pa_range, melt_viscosities_lookup, pvec(i), 'linear');
    v_vec(i) = interp1(P_Pa_range, v_lookup, pvec(i), 'linear');
    dvdp_vec(i) = interp1(P_Pa_range, dvdp_lookup, pvec(i), 'linear');

    % calculate bulk density depth vector
    if ngas_vec(i) > 0
        rhovec_bulk(i) = ((ngas_vec(i) ./ rhovec_gas(i)) + ((1.- ngas_vec(i))./rhovec_melt(i))).^(-1); 
    else 
        rhovec_bulk(i) = rhovec_melt(i);
    end


end


% apparent viscosity vector (these are all shear viscosities)
porosityvec = (rhovec_melt - rhovec_bulk) ./ (rhovec_melt - rhovec_gas); % eqn 40, CK22
muvec = muvec_melt ./ (1 - porosityvec); % eqn 40, CK22

% calculate mixture soundspeed
cvec = sqrt((gradient(rhovec_bulk)./gradient(pvec)).^(-1));


% alternative calculation for bulk soundspeed (doesn't seem very correct)
% K_vec = -(1./v_vec).*dvdp_vec; % mixture compressibility
% % bulk soundspeed
% cvec = (rhovec_bulk.^(-1).*(K_vec.*((-1)+ngas_vec).*pvec+(-1).*ngas_vec.*(K_vec+pvec).*params.R.*rhovec_bulk.* ...
%     params.T).^2.*((-1).*K_vec.*((-1)+ngas_vec).*pvec.^2+ngas_vec.*(K_vec+pvec).^2.*params.R.*rhovec_bulk.* ...
%     params.T).^(-1)).^(1/2);  


% save background state
% in Chao's model, z=0 is bottom of column, so we flip all of these
% profiles here
bg.zvec_col = zvec_col;
bg.pvec = flip(pvec);
bg.rhovec = flip(rhovec_bulk);
bg.muvec = flip(muvec);
bg.cvec = flip(cvec);
bg.rhovec_melt = flip(rhovec_melt);
bg.ngas_vec = flip(ngas_vec);
% bg.K_vec = flip(K_vec);

%% plot background state %%

% % % plot colors % % %
red = "#BF4539";
orange = "#D9863D";
green = "#BCBF65";
blue = "#51A6A6";
% % % % % % % % % % % %

figure(1)
subplot(1,4,1)
plot(bg.rhovec, bg.zvec_col, 'color', red, 'LineWidth', 1.5)
ylabel('height in column m', "FontSize", 12);
xlabel('bulk density kg/m3', "FontSize", 12);
grid
% set(gca, 'YDir','reverse')

subplot(1,4,2)
plot(bg.pvec, bg.zvec_col, 'color', orange, 'LineWidth', 1.5)
xlabel('pressure Pa', "FontSize", 12);
grid

subplot(1,4,3)
plot(bg.ngas_vec, bg.zvec_col, 'color', green, 'LineWidth', 1.5)
xlabel('n_{gas}', "FontSize", 12);
grid

subplot(1,4,4)
plot(bg.cvec, bg.zvec_col, 'color', blue, 'LineWidth', 1.5)
xlabel('soundspeed m/s', "FontSize", 12);
grid

end

%% RHS ODE, solve for pressure as a function of z %%

function RHS = RHS(t,y,params, data)

% unpack data struture
pressure_Pa_table = data.Pressure_Pa;
melt_density_table = data.rho_melt_kgm3; 
gas_density_table = data.rho_gas_kgm3; 

% extract table values
P_Pa_range = pressure_Pa_table(1:end);
gas_densities_lookup = gas_density_table(1:end);
melt_densities_lookup = melt_density_table(1:end);

P=y;

% 1D Interpolation to find densities at P
% initialize vectors
rho_gas = nan(size(P));
rho_melt = nan(size(P));

for i = 1:length(P)
    rho_gas(i) = interp1(P_Pa_range, gas_densities_lookup, P(i), 'linear', 'extrap');
    rho_melt(i) = interp1(P_Pa_range, melt_densities_lookup, P(i), 'linear', 'extrap');
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
