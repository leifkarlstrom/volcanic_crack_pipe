function rhovec_gas = DZ2006EOS(P, T, H2Ofrac)

    % Inputs:
    % P: pressure vector (Pa)
    % T: magma temperature (K)
    % H2Ofrac: H2O molar fraction
    % Outputs:
    % rhovec_gas: vector of gas density (kg/m^3)

    %YAXIS_P_MPa_range = np.linspace(0.1, 200, 100)
    %ZAXIS_T_C_range = np.linspace(1000, 1500, 100)
    %XAXIS_XH2O_range = np.linspace(0, 1, 100)

    % values based on Gas Equation of State (Duan & Zhang, 2006)
    % EOS of the H2O, CO2, and H2O–CO2 systems up to 10 GPa and 2573.15 K

        % interpolation range for current lookup table
        % P: 0.1 - 200 MPa
        % T: 1000 - 1500 deg C
        % XH2O: 0-1

    % lookup table generated from DiadFit by Penny Weiser
        % Bounds for lookup table: 
        % Pressure ()
        % Temp ()
        % H20 fraction (0 - 1)

    % Convert P (Pa) to P (MPa)
    P_MPa = P / 1E6;

    % Load DZ2006EOS lookup table from the .mat file
    data = load('D&ZEOS_lookuptable_50x50x10XH2O_MPa_kgm3_K.mat');

    % Extract the 3D array from the loaded data
    densities_array = data.densities_with_axes;

    % Extract axis values
    XH2O_range = densities_array(2:end, 1, 1);   % XH2O axis
    P_MPa_range = densities_array(1, 2:end, 1);  % Pressure axis
    T_K_range = densities_array(1, 1, 2:end);    % Temperature axis

    % Extract density values
    densities = densities_array(2:end, 2:end, 2:end);

    [P_grid, XH2O_grid, T_grid] = meshgrid(P_MPa_range, XH2O_range, T_K_range);


    % Initialize output vector
    rhovec_gas = nan(size(P_MPa));

    % Loop over the pressure values
    for i = 1:length(P_MPa)
        % Perform interpolation for each pressure
        rhovec_gas(i) = interp3(P_grid, XH2O_grid, T_grid, densities, ...
                                P_MPa(i), H2Ofrac, T, 'linear');
    end

    % Check for extrapolation warnings
    if any(isnan(rhovec_gas))
        warning('Some requested values are outside the interpolation range.');
    end
end

