function bg = magmastatic(params)

% MAGMASTATIC computes magmastatic background state
    % created by Keel Wilde, Feb. 28, 2024

%   variables
    
    % n_tot: total volatile contents (H2O + CO2)        [mass fraction]
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

% 
% %% import the PetThermoTool (alphaMELTS) python script as a module
% % needs to be done inside of every function that calls it
% run_PTT_module = py.importlib.import_module('run_PTT');
% py.importlib.reload(run_PTT_module);  % Reload in case of changes

%% 

% prepare MELTS lookup tables
data = load(sprintf('alphaMELTS_bulk_density_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
bulk_density_table = data.bulk_density_kg_m3; 
data = load(sprintf('alphaMELTS_melt_density_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
melt_density_table = data.melt_density_kg_m3;
data = load(sprintf('alphaMELTS_gas_density_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
gas_density_table = data.gas_density_kg_m3;
data = load(sprintf('alphaMELTS_melt_viscosity_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
melt_viscosity_table = data.melt_viscosity_Pas;
data = load(sprintf('alphaMELTS_n_gas_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
ngas_table = data.n_gas;
data = load(sprintf('alphaMELTS_v_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
v_table = data.volume_mixture_cm3;
data = load(sprintf('alphaMELTS_dvdp_PressurePa_H2Ofrac_%.1fC_ntot_%.3f.mat', params.T, params.n_tot));
dvdp_table = data.dvdp_mixture;


% convert Temp celcius to temp kelvin
params.T = params.T + 273.15;

% define depth vector (full column plus res)(not used here currently)
% zvec = (0:params.dz:2*params.Rres+params.Lcol);         
 
zvec_col = params.z;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(@(t,y) RHS(t,y,params, bulk_density_table),[0 params.Lcol],[params.p_atm],options);      % solve ODE

%evaluate solution and derivatives at high order
[pvec,dPdz]=deval(zvec_col,sol); 

% re-orient pvec
pvec = pvec';

% initialize depth vectors
ngas_vec = NaN(length(pvec),1);
rhovec_gas = NaN(length(pvec),1);
rhovec_melt = NaN(length(pvec),1);
rhovec_bulk = NaN(length(pvec),1);


switch params.gasexsolve

    % Henryʻs law
    case 'henrys law'

        for i = 1:length(pvec)

            % melt density
            rhovec_melt(i) = params.rho0 * (1+ pvec(i)/params.K_melt);

            if pvec(i) < (params.n_tot / params.sigma)^2
                ngas_vec(i) = (params.n_tot-params.sigma .*pvec(i).^(1/2)); 
        
        
                % depth vector of gas density
                switch params.gasEOS
        
                    % ideal gas law
                    case 'ideal'
                
                        rhovec_gas(i) = pvec(i) ./ (params.R .*params.T);
                
                    case 'H2OCO2'
                
                        rhovec_gas(i) = DZ2006EOS(pvec(i),params.T,params.H2Ofrac);
                end
             
                rhovec_bulk(i) = ((ngas_vec(i) ./ rhovec_gas(i)) + ((1.-ngas_vec(i))./rhovec_melt(i))).^(-1);
        
            else 
                % no gas exsolved
                ngas_vec(i) = 0;
                disp("No gas exsolved at " + zvec_col(i))
        
                % so no gas density
                rhovec_gas(i) = NaN;
                % % depth vector of melt density
                  
                rhovec_bulk(i) = rhovec_melt(i);
        
            end
        end

    case 'alphaMELTS'
     
        % call lookup tables
        % extract axis values (we can just do this once)
        H2Ofrac_range = bulk_density_table(2:end, 1);   % H2Ofrac axis
        P_Pa_range = bulk_density_table(1, 2:end);  % Pressure axis
        
        % extract table values
        bulk_densities_lookup = bulk_density_table(2:end, 2:end);
        melt_densities_lookup = melt_density_table(2:end, 2:end);
        gas_densities_lookup = gas_density_table(2:end, 2:end);
        melt_viscosities_lookup = melt_viscosity_table(2:end, 2:end);
        ngas_lookup = ngas_table(2:end, 2:end);
        v_lookup = v_table(2:end, 2:end);
        dvdp_lookup = dvdp_table(2:end, 2:end);
        
        
        % refine table for specific H2Ofrac value
        [~, row_idx] = min(abs(H2Ofrac_range - params.H2Ofrac));
        bulk_densities_lookup = bulk_densities_lookup(row_idx, :); 
        melt_densities_lookup = melt_densities_lookup(row_idx, :);
        gas_densities_lookup = gas_densities_lookup(row_idx, :);
        melt_viscosities_lookup = melt_viscosities_lookup(row_idx, :);
        ngas_lookup = ngas_lookup(row_idx, :);
        v_lookup = v_lookup(row_idx, :);
        dvdp_lookup = dvdp_lookup(row_idx, :);


        % interpolate lookup table to find rho for a given P
        % Initialize output vector
        rhovec_bulk = nan(size(pvec));
        rhovec_melt = nan(size(pvec));
        rhovec_gas = nan(size(pvec));
        muvec_melt = nan(size(pvec));
        ngas_vec = nan(size(pvec));
        v_vec = nan(size(pvec));
        dvdp_vec = nan(size(pvec));
    
        % loop over the pressure values
        for i = 1:length(pvec)
            % 1D Interpolation to find depth vectors
            rhovec_bulk(i) = interp1(P_Pa_range, bulk_densities_lookup, pvec(i), 'linear');
            rhovec_melt(i) = interp1(P_Pa_range, melt_densities_lookup, pvec(i), 'linear');
            rhovec_gas(i) = interp1(P_Pa_range, gas_densities_lookup, pvec(i), 'linear');
            muvec_melt(i) = interp1(P_Pa_range, melt_viscosities_lookup, pvec(i), 'linear');
            ngas_vec(i) = interp1(P_Pa_range, ngas_lookup, pvec(i), 'linear');
            v_vec(i) = interp1(P_Pa_range, v_lookup, pvec(i), 'linear');
            dvdp_vec(i) = interp1(P_Pa_range, dvdp_lookup, pvec(i), 'linear');
        
        end
      
        % check for extrapolation warnings
        if any(isnan(rhovec_bulk))
            warning('Some requested values are outside the interpolation range.');
        end
  

        % apparent viscosity vector (these are all shear viscosities)
        porosityvec = (rhovec_melt - rhovec_bulk) ./ (rhovec_melt - rhovec_gas); % eqn 40, CK22
        muvec = muvec_melt ./ (1 - porosityvec); % eqn 40, CK22
        bg.muvec = flip(muvec);
          
        



end





% this is currently not used but also needs to be re-written for the case
% of H2OCO2 
% drhog_dp = 1/(params.R .*params.T);

% % depth vector of melt density
% rhovec_melt = params.rho0.*(1+ pvec./params.K_melt);

%currently unused
% drhom_dp = params.rho0/params.K_melt;

% % depth vector of bulk density
% for i = 1:length(nvec_gas)
%     if nvec_gas(i) > 0
%         rhovec_bulk(i) = ((nvec_gas ./ rhovec_gas) + ((1.-nvec_gas)./rhovec_melt)).^(-1); 
%     
%     else 
%         rhovec_bulk = rhovec_melt;
%     end
% end

% depth vector of soundspeed
% drhodP = gradient(rhovec_bulk) ./ gradient(pvec);
% cvec = drhodP.^(-1/2);
Beta_vec = -(1./v_vec).*dvdp_vec; % mixture compressibility
K_vec = 1./Beta_vec; %bulk mod (Pa)

cvec = (rhovec_bulk.^(-1).*(K_vec.*((-1)+ngas_vec).*pvec+(-1).*ngas_vec.*(K_vec+pvec).*params.R.*rhovec_bulk.* ...
    params.T).^2.*((-1).*K_vec.*((-1)+ngas_vec).*pvec.^2+ngas_vec.*(K_vec+pvec).^2.*params.R.*rhovec_bulk.* ...
    params.T).^(-1)).^(1/2);  

% save background state
% in Chao's model, z=0 is bottom of column, so we flip all of these
% profiles here
bg.zvec_col = zvec_col;
bg.pvec = flip(pvec);
bg.rhovec = flip(rhovec_bulk);
bg.cvec = flip(cvec);
bg.rhovec_melt = flip(rhovec_melt);
bg.ngas_vec = flip(ngas_vec);
bg.K_vec = flip(K_vec);


% figure(1)
% subplot(1,4,1)
% plot(bg.rhovec, bg.zvec_col)
% ylabel('depth in column (m)', "FontSize", 12);
% xlabel('bulk density (kg/m3 ?)', "FontSize", 12);
% % set(gca, 'YDir','reverse')
% 
% 
% subplot(1,4,2)
% plot(bg.pvec, bg.zvec_col)
% ylabel('depth in column (m)', "FontSize", 12);
% xlabel('pressure', "FontSize", 12);
% % set(gca, 'YDir','reverse')
% 
% subplot(1,4,3)
% plot(bg.ngas_vec, bg.zvec_col)
% ylabel('depth in column (m)', "FontSize", 12);
% xlabel('ngas', "FontSize", 12);
% % set(gca, 'YDir','reverse')
% 
% subplot(1,4,4)
% plot(bg.cvec, bg.zvec_col)
% ylabel('depth in column (m)', "FontSize", 12);
% xlabel('cvec', "FontSize", 12);
% set(gca, 'YDir','reverse')

end

%% RHS ODE, solve for pressure as a function of z %%
function RHS = RHS(t,y,params, bulk_density_table)

% call lookup table
% extract axis values
H2Ofrac_range = bulk_density_table(2:end, 1);   % H2Ofrac axis
P_Pa_range = bulk_density_table(1, 2:end);  % Pressure axis
% extract table values
bulk_densities_lookup = bulk_density_table(2:end, 2:end);

% refine table for specific H2Ofrac value
[~, row_idx] = min(abs(H2Ofrac_range - params.H2Ofrac));
bulk_densities_lookup = bulk_densities_lookup(row_idx, :); 



P=y;

% if you use this you need to update melt bulk mod
switch params.gasexsolve

    % Henryʻs law
    case 'henrys law'

        if P(1) < (params.n_tot / params.sigma)^2
            n_gas = (params.n_tot-params.sigma .* P.^(1/2)); 

            % EOS for gas density
            switch params.gasEOS

                % ideal gas law
                case 'ideal'

                    rho_gas = P ./ (params.R .*params.T);

                case 'H2OCO2'

                    rho_gas = DZ2006EOS(P,params.T,params.H2Ofrac);
            end

        else 
             % no gas exsolved
            n_gas = 0;

            % so no gas density
            rho_gas = NaN;



        end

        % melt density
        rho_melt = params.rho0.*(1+ P./params.K_melt);

        % rho_melt = PTT_results.melt_density;

        %bulk density
        if n_gas > 0
            rho = ((n_gas ./ rho_gas) + ((1.-n_gas)./rho_melt)).^(-1); 

        else 
            rho = rho_melt;

        end


    case 'alphaMELTS'
        % interpolate lookup table to find rho for a given P
        % Initialize output vector
        rho = nan(size(P));

        % loop over the pressure values
        for i = 1:length(P)
            % 1D Interpolation to find bulk density at P
            rho(i) = interp1(P_Pa_range, bulk_densities_lookup, P(i), 'linear');

        end
       % check for extrapolation warnings
        if any(isnan(rho))
            warning('Some requested values are outside the interpolation range.');
        end


        



 end

RHS = (rho.*params.g );


end
