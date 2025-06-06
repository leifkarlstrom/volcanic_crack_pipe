function mu_melt = melt_viscosity(params)


%% composition
wtp_sum_nonH2O = 50.6+2.3+12.9+10.56+0.87+0.19+8.95+10.6+2.19+0.38+0.23+0;
wtp_scale = (100)/wtp_sum_nonH2O; %adjust so that % adds to 100
wtp_SiO2 = wtp_scale*50.6;
wtp_TiO2 = wtp_scale*2.3;
wtp_Al2O3 = wtp_scale*12.9;
wtp_FeO = wtp_scale*10.56; 
wtp_Fe2O3 = wtp_scale*0.87;
% wtp_FeO_Total = wtp_FeO + wtp_Fe2O3;
wtp_MnO = wtp_scale*0.19;
wtp_MgO = wtp_scale*8.95;
wtp_CaO = wtp_scale*10.6;
wtp_Na2O = wtp_scale*2.19;
wtp_K2O = wtp_scale*0.38;
wtp_P2O5 = wtp_scale*0.23;
wtp_F2O_1 = wtp_scale*0; %assuming

%% molar weights (kg/mol)
M_SiO2 = 60.09/1000;
M_TiO2 = 79.88/1000;
M_Al2O3 = 101.96/1000;
M_Fe2O3 = 159.7/1000;
M_FeO = 71.85/1000;
M_MnO = 70.94/1000;
M_MgO = 40.31/1000;
M_CaO = 56.08/1000;
M_Na2O = 61.98/1000;
M_K2O = 94.2/1000;
M_P2O5 = 141.94/1000;
M_H2O = 18.02/1000;
M_F2O_1 = 54.0/1000;
M_CO2 = 44.01/1000;

%% mol/kg
mol_kg_SiO2 = wtp_SiO2/100/(M_SiO2);
mol_kg_TiO2 = wtp_TiO2/100/(M_TiO2);
mol_kg_Al2O3 = wtp_Al2O3/100/(M_Al2O3);
mol_kg_Fe2O3 = wtp_Fe2O3/100/(M_Fe2O3);
mol_kg_FeO = wtp_FeO/100/(M_FeO);
%mol_kg_FeO_Total = mol_kg_Fe2O3 + mol_kg_FeO;
mol_kg_MnO = wtp_MnO/100/(M_MnO);
mol_kg_MgO = wtp_MgO/100/(M_MgO);
mol_kg_CaO = wtp_CaO/100/(M_CaO);
mol_kg_Na2O = wtp_Na2O/100/(M_Na2O);
mol_kg_K2O = wtp_K2O/100/(M_K2O);
mol_kg_P2O5 = wtp_P2O5/100/(M_P2O5);
% mol_kg_H2O = wtp_H2O/100/(M_H2O/1000);
mol_kg_F2O_1 = wtp_F2O_1/100/(M_F2O_1);

mol_kg_total_dry = mol_kg_SiO2 + mol_kg_TiO2 + mol_kg_Al2O3 + mol_kg_FeO + mol_kg_Fe2O3 +...
    +mol_kg_MnO + mol_kg_MgO + mol_kg_CaO + mol_kg_Na2O + mol_kg_K2O + mol_kg_P2O5;

%% convert ppm to wt%
% H2Osolwtp = H2Osolppm/1E6*M_H2O*mol_kg_total*100;
CO2solwtp_table = CO2solppm_table/1E6*M_CO2*mol_kg_total_dry*100;

CO2_condtop_wtp = CO2_condtop_ppm/1E6*M_CO2*mol_kg_total_dry*100;
CO2_condbot_wtp = CO2_condbot_ppm/1E6*M_CO2*mol_kg_total_dry*100;
CO2_laketop_wtp = CO2_laketop_ppm/1E6*M_CO2*mol_kg_total_dry*100;
CO2_lakebot_wtp = CO2_lakebot_ppm/1E6*M_CO2*mol_kg_total_dry*100;

%% depth vectors
zvec = (0:dz:2*Rres+Hcolumn+dz)';
Pvec = NaN(size(zvec));
H2O_gas_vec = NaN(size(zvec));
CO2_gas_vec = NaN(size(zvec));

H2O_gas_vec_p = NaN(size(zvec));
CO2_gas_vec_p = NaN(size(zvec));

gas_massfrac_vec = NaN(size(zvec));
H2O_gas_mol_vec = NaN(size(zvec));
CO2_gas_mol_vec = NaN(size(zvec)); 
H2O_gas_molfrac_vec = NaN(size(zvec));

gas_vol_vec = NaN(size(zvec));
melt_vol_vec = NaN(size(zvec));
porosity_vec = NaN(size(zvec));
density_vec = NaN(size(zvec));
gas_density_vec = NaN(size(zvec));
dwtp_H2O_gas_dP_vec = NaN(size(zvec));

dwtp_CO2_gas_dP_vec = NaN(size(zvec));
dgasmassfrac_dP_vec = NaN(size(zvec));
dgasdensity_dP_vec = NaN(size(zvec));
ddensity_dP_vec = NaN(size(zvec));

wtp_H2O_dis_vec = NaN(size(zvec));
mol_kg_H2O_vec = NaN(size(zvec));
mol_total_vec = NaN(size(zvec));
SiO2_vec = NaN(size(zvec));
TiO2_vec = NaN(size(zvec));
Al2O3_vec = NaN(size(zvec));
Fe2O3_vec = NaN(size(zvec));
FeO_vec = NaN(size(zvec));
FeO_Total_vec = NaN(size(zvec));
MnO_vec = NaN(size(zvec));
MgO_vec = NaN(size(zvec));
CaO_vec = NaN(size(zvec));
Na2O_vec = NaN(size(zvec));
K2O_vec = NaN(size(zvec));
P2O5_vec = NaN(size(zvec));
H2Omol_vec = NaN(size(zvec));
F2O_1_vec = NaN(size(zvec));
melt_kg_mol_vec = NaN(size(zvec));
melt_density_vec = NaN(size(zvec));

%% H2O-CO2
CO2_vec = CO2_condbot_wtp+zeros(size(zvec));
CO2_vec(zvec<Hcolumn & zvec>=Hlake) = CO2_condtop_wtp +...
    (CO2_condbot_wtp - CO2_condtop_wtp)/(Hcolumn-Hlake)*...
    (zvec(zvec<Hcolumn & zvec>=Hlake) - Hlake);
CO2_vec(zvec<Hlake) = CO2_laketop_wtp +...
    (CO2_lakebot_wtp - CO2_laketop_wtp)/(Hlake)*...
    zvec(zvec<Hlake);
H2O_vec = H2O_condbot_wtp+zeros(size(zvec));
H2O_vec(zvec<Hcolumn & zvec>=Hlake) = H2O_condtop_wtp +...
    (H2O_condbot_wtp - H2O_condtop_wtp)/(Hcolumn-Hlake)*...
    (zvec(zvec<Hcolumn & zvec>=Hlake) - Hlake);
H2O_vec(zvec<Hlake) = H2O_laketop_wtp +...
    (H2O_lakebot_wtp - H2O_laketop_wtp)/(Hlake)*...
    zvec(zvec<Hlake);
TK_vec = TK_condbot+zeros(size(zvec));
TK_vec(zvec<Hcolumn & zvec>=Hlake) = TK_condtop +...
    (TK_condbot - TK_condtop)/(Hcolumn-Hlake)*...
    (zvec(zvec<Hcolumn & zvec>=Hlake) - Hlake);
TK_vec(zvec<Hlake) = TK_laketop +...
    (TK_lakebot - TK_laketop)/(Hlake)*...
    zvec(zvec<Hlake);
for i = 1:length(zvec)
    if i == 1
        Pvec(i) = 1E5; %atmospheric pressure
    else
        Pvec(i) = Pvec(i-1) + dz*density_vec(i-1)*g;
    end
    
    %% find wtp exsolved gas at each pressure
    %mol gas/kg melt
%     H2Osolwtp_temp = interp2(H2O_gas_molfrac_table, P_table, ...
%         H2Osolwtp_table, H2O_gas_molfrac_table, Pvec(i));
%     CO2solwtp_temp = interp2(H2O_gas_molfrac_table, P_table, ...
%         CO2solwtp_table, H2O_gas_molfrac_table, Pvec(i));
    Pin_high = find(P_table>Pvec(i), 1);
    if Pin_high == 1
        Pin_high = 2;
    elseif isempty(Pin_high)
        error('no pressure found in table')
    end
    stepfrac = (Pvec(i)-P_table(Pin_high-1))/(P_table(Pin_high)-P_table(Pin_high-1));
    H2Osolwtp_temp = H2Osolwtp_table(Pin_high-1,:) + ...
        stepfrac*(H2Osolwtp_table(Pin_high,:)-H2Osolwtp_table(Pin_high-1,:));
    CO2solwtp_temp = CO2solwtp_table(Pin_high-1,:) + ...
        stepfrac*(CO2solwtp_table(Pin_high,:)-CO2solwtp_table(Pin_high-1,:));      
    temp_H2O_gas_mol = (H2O_vec(i)-H2Osolwtp_temp)/M_H2O;
    temp_CO2_gas_mol = (CO2_vec(i)-CO2solwtp_temp)/M_CO2;
    some_H2O_gas = temp_H2O_gas_mol > 0;
    some_CO2_gas = temp_CO2_gas_mol > 0;
    temp_H2O_gas_molfrac = NaN(size(temp_H2O_gas_mol));
    temp_H2O_gas_molfrac(some_H2O_gas & some_CO2_gas)...
        = temp_H2O_gas_mol(some_H2O_gas & some_CO2_gas)...
        ./(temp_CO2_gas_mol(some_H2O_gas & some_CO2_gas)...
        + temp_H2O_gas_mol(some_H2O_gas & some_CO2_gas));
    temp_H2O_gas_molfrac(some_H2O_gas & ~some_CO2_gas) = 1;
    temp_H2O_gas_molfrac(~some_H2O_gas & some_CO2_gas) = 0;

    %find first H2O gas ratio bigger than actual
    bigind = find(temp_H2O_gas_molfrac >= H2O_gas_molfrac_table & ~isnan(temp_H2O_gas_molfrac),1);
    %find first H2O gas ratio smaller than actual
    smallind = find(temp_H2O_gas_molfrac <= H2O_gas_molfrac_table & ~isnan(temp_H2O_gas_molfrac),1);
    %take whichever index appropriate
    if isempty(bigind) && isempty(smallind)
        ind = NaN; %no valid points
    elseif isempty(bigind)
        if smallind ~= 1
            ind = NaN;
        end
    elseif isempty(smallind)
        if bigind ~= length(temp_H2O_gas_molfrac)
            ind = NaN;
        end
    else
        ind = max([bigind,smallind]);
        if ind == 1
            ind = 2;%so don't have to add more cases to interpolation below
        end
    end


    if ~isnan(ind)
        if ~isnan(temp_H2O_gas_molfrac(ind-1)) && ~isnan(temp_H2O_gas_molfrac(ind))
            %linear interp to optimal gasmolfrac
            target_H2O_gas_molfrac = H2O_gas_molfrac_table(ind-1)...
                + (temp_H2O_gas_molfrac(ind-1) - H2O_gas_molfrac_table(ind-1))...
                /(1 - (temp_H2O_gas_molfrac(ind)-temp_H2O_gas_molfrac(ind-1))...
                /delta_H2O_gas_molfrac_table);
            %how far to step from undershoot value to overshoot value
            molfrac_step = (target_H2O_gas_molfrac - H2O_gas_molfrac_table(ind-1));
            %linear interp exsolved gasses at target gasmolfrac
            stepfrac = molfrac_step/delta_H2O_gas_molfrac_table;
            if stepfrac > 1
                %'error'
                stepfrac = 1;
            elseif stepfrac < 0
                %'error'
                stepfrac = 0;
            end
            H2O_gas_vec(i) = H2O_vec(i) - (H2Osolwtp_temp(ind-1) + stepfrac...
                *(H2Osolwtp_temp(ind)-H2Osolwtp_temp(ind-1)));
            CO2_gas_vec(i) = CO2_vec(i) - (CO2solwtp_temp(ind-1) + stepfrac...
                *(CO2solwtp_temp(ind)-CO2solwtp_temp(ind-1)));                           
        else
            %just use overshoot value
            H2O_gas_vec(i) = H2O_vec(i) - H2Osolwtp_temp(ind);
            CO2_gas_vec(i) = CO2_vec(i) - CO2solwtp_temp(ind);
        end
        %seems to be neccessary due to interp inaccuracy
        if H2O_gas_vec(i) < 0
            H2O_gas_vec(i) = 0;
        end
        if CO2_gas_vec(i) < 0
            CO2_gas_vec(i) = 0;
        end
    else
        H2O_gas_vec(i) = 0;
        CO2_gas_vec(i) = 0;
    end    
    %necessary due to floating point or interp error
    if H2O_gas_vec(i) > H2O_vec(i)
        H2O_gas_vec(i) = H2O_vec(i);
    end
    if CO2_gas_vec(i) > CO2_vec(i)
        CO2_gas_vec(i) = CO2_vec(i);
    end


    %% find wtp exsolved gas at each [pressure + delta P]  
    %mol gas/kg melt
%     H2Osolwtp_temp = interp2(H2O_gas_molfrac_table, P_table, ...
%         H2Osolwtp_table, H2O_gas_molfrac_table, Pvec(i)+deltaP);
%     CO2solwtp_temp = interp2(H2O_gas_molfrac_table, P_table, ...
%         CO2solwtp_table, H2O_gas_molfrac_table, Pvec(i)+deltaP);
    Pin_high = find(P_table>Pvec(i)+deltaP, 1);
    if Pin_high == 1
        Pin_high = 2;
    elseif isempty(Pin_high)
        error('no pressure found in table')
    end
    stepfrac = (Pvec(i)+deltaP-P_table(Pin_high-1))/(P_table(Pin_high)-P_table(Pin_high-1));
    H2Osolwtp_temp = H2Osolwtp_table(Pin_high-1,:) + ...
        stepfrac*(H2Osolwtp_table(Pin_high,:)-H2Osolwtp_table(Pin_high-1,:));
    CO2solwtp_temp = CO2solwtp_table(Pin_high-1,:) + ...
        stepfrac*(CO2solwtp_table(Pin_high,:)-CO2solwtp_table(Pin_high-1,:)); 
    temp_H2O_gas_mol = (H2O_vec(i)-H2Osolwtp_temp)/M_H2O;
    temp_CO2_gas_mol = (CO2_vec(i)-CO2solwtp_temp)/M_CO2;
    some_H2O_gas = temp_H2O_gas_mol > 0;
    some_CO2_gas = temp_CO2_gas_mol > 0;
    temp_H2O_gas_molfrac = NaN(size(temp_H2O_gas_mol));
    temp_H2O_gas_molfrac(some_H2O_gas & some_CO2_gas)...
        = temp_H2O_gas_mol(some_H2O_gas & some_CO2_gas)...
        ./(temp_CO2_gas_mol(some_H2O_gas & some_CO2_gas)...
        + temp_H2O_gas_mol(some_H2O_gas & some_CO2_gas));
    temp_H2O_gas_molfrac(some_H2O_gas & ~some_CO2_gas) = 1;
    temp_H2O_gas_molfrac(~some_H2O_gas & some_CO2_gas) = 0;

    %find first H2O gas ratio bigger than actual
    bigind = find(temp_H2O_gas_molfrac >= H2O_gas_molfrac_table & ~isnan(temp_H2O_gas_molfrac),1);
    %find first H2O gas ratio smaller than actual
    smallind = find(temp_H2O_gas_molfrac <= H2O_gas_molfrac_table & ~isnan(temp_H2O_gas_molfrac),1);
    %take whichever index appropriate
    if isempty(bigind) && isempty(smallind)
        ind = NaN; %no valid points
    elseif isempty(bigind)
        if smallind ~= 1
            ind = NaN;
        end
    elseif isempty(smallind)
        if bigind ~= length(temp_H2O_gas_molfrac)
            ind = NaN;
        end
    else
        ind = max([bigind,smallind]);
        if ind == 1
            ind = 2;%so don't have to add more cases to interpolation below
        end
    end


    if ~isnan(ind)
        if ~isnan(temp_H2O_gas_molfrac(ind-1)) && ~isnan(temp_H2O_gas_molfrac(ind))
            %linear interp to optimal gasmolfrac
            target_H2O_gas_molfrac = H2O_gas_molfrac_table(ind-1)...
                + (temp_H2O_gas_molfrac(ind-1) - H2O_gas_molfrac_table(ind-1))...
                /(1 - (temp_H2O_gas_molfrac(ind)-temp_H2O_gas_molfrac(ind-1))...
                /delta_H2O_gas_molfrac_table);
            %how far to step from undershoot value to overshoot value
            molfrac_step = (target_H2O_gas_molfrac - H2O_gas_molfrac_table(ind-1));
            %linear interp exsolved gasses at target gasmolfrac
            stepfrac = molfrac_step/delta_H2O_gas_molfrac_table;
            if stepfrac > 1
                %'error'
                stepfrac = 1;
            elseif stepfrac < 0
                %'error'
                stepfrac = 0;
            end
            H2O_gas_vec_p(i) = H2O_vec(i) - (H2Osolwtp_temp(ind-1) + stepfrac...
                *(H2Osolwtp_temp(ind)-H2Osolwtp_temp(ind-1)));
            CO2_gas_vec_p(i) = CO2_vec(i) - (CO2solwtp_temp(ind-1) + stepfrac...
                *(CO2solwtp_temp(ind)-CO2solwtp_temp(ind-1)));                           
        else
            %just use overshoot value
            H2O_gas_vec_p(i) = H2O_vec(i) - H2Osolwtp_temp(ind);
            CO2_gas_vec_p(i) = CO2_vec(i) - CO2solwtp_temp(ind);
        end
        %seems to be neccessary due to interp inaccuracy
        if H2O_gas_vec_p(i) < 0
            H2O_gas_vec_p(i) = 0;
        end
        if CO2_gas_vec_p(i) < 0
            CO2_gas_vec_p(i) = 0;
        end
    else
        H2O_gas_vec_p(i) = 0;
        CO2_gas_vec_p(i) = 0;
    end  
    %necessary due to floating point or interp error
    if H2O_gas_vec_p(i) > H2O_vec(i)
        H2O_gas_vec_p(i) = H2O_vec(i);
    end
    if CO2_gas_vec_p(i) > CO2_vec(i)
        CO2_gas_vec_p(i) = CO2_vec(i);
    end
    
    
    %% melt density
    wtp_H2O_dis_vec(i) = H2O_vec(i) - H2O_gas_vec(i);
    mol_kg_H2O_vec(i) = wtp_H2O_dis_vec(i)/100/(M_H2O);
    %mol% (needed for melt visc model)
    mol_total_vec(i) = mol_kg_SiO2 + mol_kg_TiO2 + mol_kg_Al2O3 + mol_kg_FeO + mol_kg_Fe2O3 +...
        +mol_kg_MnO + mol_kg_MgO + mol_kg_CaO + mol_kg_Na2O + mol_kg_K2O + mol_kg_P2O5 + mol_kg_H2O_vec(i);
    SiO2_vec(i) = 100*mol_kg_SiO2/mol_total_vec(i);
    TiO2_vec(i) = 100*mol_kg_TiO2/mol_total_vec(i);
    Al2O3_vec(i) = 100*mol_kg_Al2O3/mol_total_vec(i);
    Fe2O3_vec(i) = 100*mol_kg_Fe2O3/mol_total_vec(i);
    FeO_vec(i) = 100*mol_kg_FeO/mol_total_vec(i);
    FeO_Total_vec(i) = Fe2O3_vec(i) + FeO_vec(i);
    MnO_vec(i) = 100*mol_kg_MnO/mol_total_vec(i);
    MgO_vec(i) = 100*mol_kg_MgO/mol_total_vec(i);
    CaO_vec(i) = 100*mol_kg_CaO/mol_total_vec(i);
    Na2O_vec(i) = 100*mol_kg_Na2O/mol_total_vec(i);
    K2O_vec(i) = 100*mol_kg_K2O/mol_total_vec(i);
    P2O5_vec(i) = 100*mol_kg_P2O5/mol_total_vec(i);
    H2Omol_vec(i) = 100*mol_kg_H2O_vec(i)/mol_total_vec(i);
    F2O_1_vec(i) = 100*mol_kg_F2O_1/mol_total_vec(i);
    
    %%%%melt density model (Lange and Carmichael 1987)
    dVdP_SiO2 = -1.73E-15;
    dVdP_TiO2 = -2.20E-15;
    dVdP_Al2O3 = -1.82E-15;
    dVdP_FeO = -0.54E-15;
    dVdP_MgO = 0.20E-15;
    dVdP_CaO = 0.08E-15;
    dVdP_Na2O = -2.59E-15;
    dVdP_K2O = -8.79E-15;
    dVdP_Fe2O3 = 0; %assuming
    dVdP_MnO = 0; %assuming
    dVdP_P2O5 = 0;%assuming
    dVdP_H20 = -3.15E-15; %Ochs and Lange 1999

    dVdT_SiO2 = 0.00E-9;
    dVdT_TiO2 = 7.24E-9;
    dVdT_Al2O3 = 2.62E-9;
    dVdT_Fe2O3 = 9.09E-9;
    dVdT_FeO = 2.92E-9;
    dVdT_MgO = 2.62E-9;
    dVdT_CaO = 2.92E-9;
    dVdT_Na2O = 7.41E-9;
    dVdT_K2O = 5.25E-9;
    dVdT_MnO = 0; %assuming
    dVdT_P2O5 = (69-88E-6)*10^-9;%from Knoche et al 1995
    dVdT_H20 = 9.46E-9; %Ochs and Lange 1999

    %molar volume (m^3/mol)
    V_SiO2 = 26.92E-6;
    V_TiO2 = 23.89E-6;
    V_Al2O3 = 37.37E-6;
    V_Fe2O3 = 42.97E-6;
    V_FeO = 13.97E-6;
    V_MgO = 11.73E-6;
    V_CaO = 16.85E-6;
    V_Na2O = 29.51E-6;
    V_K2O = 47.07E-6;
    V_excess_Na2O_TiO2 = 20.21E-6;
    V_MnO = 13.97E-6; %assuming similar to FeO

    V_MnO = 13.97E-6; %assuming similar to FeO
    V_P2O5 = 69E-6;%from Knoche et al 1995
    V_H20 = 22.89E-6 + 500*dVdT_H20;%Ochs and Lange 1999, corrected to 1773 K

    %melt volume/mol at 1773 K and 10^5 Pa
    V_melt_1773_10_5 = 1/100*(SiO2_vec(i)*V_SiO2 + TiO2_vec(i)*V_TiO2 + Al2O3_vec(i)*V_Al2O3 + FeO_vec(i)*V_FeO +...
        Fe2O3_vec(i)*V_Fe2O3 + MnO_vec(i)*V_MnO + MgO_vec(i)*V_MgO + CaO_vec(i)*V_CaO + Na2O_vec(i)*V_Na2O +...
        K2O_vec(i)*V_K2O + P2O5_vec(i)*V_P2O5 + H2Omol_vec(i)*V_H20 + Na2O_vec(i)*TiO2_vec(i)/100*V_excess_Na2O_TiO2);

    dV_melt_dT = 1/100*...
        (SiO2_vec(i)*dVdT_SiO2 + TiO2_vec(i)*dVdT_TiO2 + Al2O3_vec(i)*dVdT_Al2O3 + FeO_vec(i)*dVdT_FeO +...
        Fe2O3_vec(i)*dVdT_Fe2O3 + MnO_vec(i)*dVdT_MnO + MgO_vec(i)*dVdT_MgO + CaO_vec(i)*dVdT_CaO + Na2O_vec(i)*dVdT_Na2O +...
        K2O_vec(i)*dVdT_K2O + P2O5_vec(i)*dVdT_P2O5 + H2Omol_vec(i)*dVdT_H20);

    dV_melt_dP = 1/100*...
        (SiO2_vec(i)*dVdP_SiO2 + TiO2_vec(i)*dVdP_TiO2 + Al2O3_vec(i)*dVdP_Al2O3 + FeO_vec(i)*dVdP_FeO +...
        Fe2O3_vec(i)*dVdP_Fe2O3 + MnO_vec(i)*dVdP_MnO + MgO_vec(i)*dVdP_MgO + CaO_vec(i)*dVdP_CaO + Na2O_vec(i)*dVdP_Na2O +...
        K2O_vec(i)*dVdP_K2O + P2O5_vec(i)*dVdP_P2O5 + H2Omol_vec(i)*dVdP_H20);

    %neglecting d^2V/dPDT
    V_melt = V_melt_1773_10_5 + (TK_vec(i) - 1773)*dV_melt_dT + (Pvec(i) - 10^5)*dV_melt_dP;

    %total kg/mol melt
    melt_kg_mol_vec(i) = 1/100*(SiO2_vec(i)*M_SiO2 + TiO2_vec(i)*M_TiO2 + Al2O3_vec(i)*M_Al2O3 + FeO_vec(i)*M_FeO +...
        Fe2O3_vec(i)*M_Fe2O3 + MnO_vec(i)*M_MnO + MgO_vec(i)*M_MgO + CaO_vec(i)*M_CaO + Na2O_vec(i)*M_Na2O +...
        K2O_vec(i)*M_K2O + P2O5_vec(i)*M_P2O5 + H2Omol_vec(i)*M_H2O);

    melt_density_vec(i) = melt_kg_mol_vec(i)/V_melt; 
    dmeltdensity_dP = -melt_kg_mol_vec(i)/V_melt^2*dV_melt_dP; %from drho/dp = d(M/V)/dP 
    
    
   
%% viscosity
%%%%melt viscosity model, Giordano et al 2008
%model B coefficients and mol% formulas
M1 = SiO2_vec + TiO2_vec; b1 = 159.6; 
M2 = Al2O3_vec; b2 = -173.3;
M3 = FeO_Total_vec + MnO_vec + P2O5_vec; b3 = 72.1;
M4 = MgO_vec; b4 = 75.7; 
M5 = CaO_vec; b5 = -39.0 ;
M6 = Na2O_vec + H2Omol_vec + F2O_1_vec; b6 = -84.1;
M7 = H2Omol_vec + F2O_1_vec + log(1 + H2Omol_vec); b7 = 141.5;
M11 = (SiO2_vec + TiO2_vec).*(FeO_Total_vec + MnO_vec + MgO_vec); b11 = -2.43;
M12 = (SiO2_vec + TiO2_vec + Al2O3_vec + P2O5_vec).*(Na2O_vec + K2O_vec + H2Omol_vec); b12 = -0.91;
M13 = (Al2O3_vec).*(Na2O_vec + K2O_vec); b13 = 17.6;

%model C coefficients and mol% formulas
N1 = SiO2_vec; c1 = 2.75;
N2 = TiO2_vec + Al2O3_vec; c2 = 15.7;
N3 = FeO_Total_vec + MnO_vec + MgO_vec; c3 = 8.3;
N4 = CaO_vec; c4 = 10.2;
N5 = Na2O_vec + K2O_vec; c5 = -12.3;
N6 = log(1 + H2Omol_vec + F2O_1_vec); c6 = -99.5;
N11 = (Al2O3_vec + FeO_Total_vec + MnO_vec + MgO_vec + CaO_vec - P2O5_vec).*(Na2O_vec + K2O_vec + H2Omol_vec + F2O_1_vec); c11 = 0.30;

%melt viscosity model
A = -4.55;
B = b1*M1 + b2*M2 + b3*M3 + b4*M4 + b5*M5 + b6*M6 + b7*M7 + b11*M11 + b12*M12 + b13*M13;
C = c1*N1 + c2*N2 + c3*N3 + c4*N4 + c5*N5 + c6*N6 + c11*N11;
mu_melt_vec = 10.^(A + B./(TK_vec - C));


















mu_melt = a;
end



