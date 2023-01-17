function output_struct = magma_prop_from_volatiles_variableA_noextra(conddip,...
    H2O_condtop_wtp,H2O_condbot_wtp,CO2_condtop_ppm,CO2_condbot_ppm,...
    H2O_laketop_wtp,H2O_lakebot_wtp,CO2_laketop_ppm,CO2_lakebot_ppm,...
    Hlake,Hcolumn,Rres,Rcondtop,Rcondbot,Rlake,...
    TK_condbot,TK_condtop,TK_lakebot,TK_laketop,...
    plotprofiles, CO2solppm_table, H2Osolwtp_table, dz)

%Inputs:
    %condip: conduit dip (deg)
    %[...] total H2O (wt.%) and CO2 (ppm) at the top and bottom of the conduit and lake
    %Hlake: lake height (m)
    %Hcolumn: total magma column height [including lake] (m)
    %Rres: reservoir radius (m)
    %[Rcondtop,Rcondbot]: conduit [top,bottom] radius (m)
    %Rlake: lake top radius (m) [currently cylindrical or funnel lake geometry]
    %omega: angular frequency (rad/s)
    %TK: Temp (K)
    %plotprofiles: true to plot magma profiles
    %[CO2solppm_table, H2Osolwtp_table] optional, can be more efficient to
        %pass these in than load each time
    %dz: depth increment (m)
    
if Hlake < 0
    Hlake = 0;
end

%% parameters
CO2_noexsol = false; %no exsolution of CO2 in compressibility calc 
%(due to low diffusivity), doesn't seem to matter in practice
Rc = 8.314;%ideal gas const
g = 9.81;
% Ptop = linspace(1,100,100)*1E5;%top reservoir pressure
%linear total H2O-CO2 between top and bottom values

% if ~exist('dz','var')
%     dz = 1; %depth step
% end
deltaP = 1E7;%pressure step for compressibility, 
% if too small will cause issues with ddensity_dP being negative at depth due to
% limited resolution if solubility tables, reccomend > 1E7

stable_error_tol = 1E0;%density stability error tolerance

%% solubility tables:
if ~exist('CO2solppm_table','var')
    %http://calcul-isto.cnrs-orleans.fr/thermodynamics/applications/#/script/js_compute~H2O-CO2%20systems~H2O%20and%20CO2%20solubility.js
    CO2solppm_table = readmatrix('new_Values_CO2solubility_ppm_comb.csv');
    H2Osolwtp_table = readmatrix('new_Values_H2Osolubility_wtp_comb.csv');
end
P_table = CO2solppm_table(2:end,1)*1E5;
H2O_gas_molfrac_table = CO2solppm_table(1,2:end);
delta_H2O_gas_molfrac_table = (H2O_gas_molfrac_table(2)-H2O_gas_molfrac_table(1));
CO2solppm_table = CO2solppm_table(2:end,2:end);
H2Osolwtp_table = H2Osolwtp_table(2:end,2:end);

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
    
    
    %% find porosity and magma density at each pressure
    wtp_H2O_dis = H2O_vec(i) - H2O_gas_vec(i);
    mol_kg_H2O_vec = wtp_H2O_dis/100/(M_H2O);
    mol_kg_total = mol_kg_SiO2 + mol_kg_TiO2 + mol_kg_Al2O3 + mol_kg_FeO + mol_kg_Fe2O3 +...
        mol_kg_MnO + mol_kg_MgO + mol_kg_CaO + mol_kg_Na2O + mol_kg_K2O + mol_kg_P2O5 + mol_kg_H2O_vec;
    gas_massfrac_vec(i) = (H2O_gas_vec(i) + CO2_gas_vec(i))/100;
    H2O_gas_mol_vec(i) = H2O_gas_vec(i)/(M_H2O*mol_kg_total*100);
    CO2_gas_mol_vec(i) = CO2_gas_vec(i)/(M_CO2*mol_kg_total*100);  
    H2O_gas_molfrac_vec(i) = H2O_gas_mol_vec(i)/(H2O_gas_mol_vec(i) + CO2_gas_mol_vec(i));
    %ideal gas law
    gas_vol_vec(i) = (H2O_gas_mol_vec(i) + CO2_gas_mol_vec(i))*Rc*TK_vec(i)/Pvec(i);
    melt_vol_vec(i) = (1-gas_massfrac_vec(i))/melt_density_vec(i);
    porosity_vec(i) = gas_vol_vec(i)/(gas_vol_vec(i) + melt_vol_vec(i));
    density_vec(i) = 1/(gas_vol_vec(i) + melt_vol_vec(i));
    %                 gas_density_vec = gas_massfrac_vec./gas_vol_vec;
    gas_density_vec(i) = (Pvec(i)*(M_H2O*H2O_gas_molfrac_vec(i)...
        + M_CO2*(1 - H2O_gas_molfrac_vec(i))))/(Rc*TK_vec(i));
    %                 %just to prevent divide-by-0 errors in compressibility calc
    if isnan(gas_density_vec(i))
        gas_density_vec(i) = Pvec(i)*M_CO2/(Rc*TK_vec(i));   
    end

    %exsolution at equilibrium solubility
    dwtp_H2O_gas_dP_vec(i) = (H2O_gas_vec_p(i) - H2O_gas_vec(i))/deltaP;
    if CO2_noexsol
        %!!!!approximate, in reality H2O also effected if no
        %CO2 exsolution, and also wouldn't be zero just less
        %than equilibrium!!!!
        dwtp_CO2_gas_dP_vec(i) = 0;
    else
        dwtp_CO2_gas_dP_vec(i) = (CO2_gas_vec_p(i) - CO2_gas_vec(i))/deltaP;
    end
    dgasmassfrac_dP_vec(i) = (dwtp_H2O_gas_dP_vec(i) + dwtp_CO2_gas_dP_vec(i))/100;
    dgasdensity_dP_vec(i) = (H2O_gas_molfrac_vec(i)*M_H2O ...
        + (1 - H2O_gas_molfrac_vec(i))*M_CO2)/(Rc*TK_vec(i)); %ideal gas law
    %just to prevent divide-by-0 errors in compressibility calc
    if isnan(dgasdensity_dP_vec(i))
        dgasdensity_dP_vec(i) = M_CO2/(Rc*TK_vec(i));
    end
    %general compressibility (derivative of eq in Kieffer 1977)
    %                 ddensity_dP_vec = (melt_density*(gas_density_vec.*dgasmassfrac_dP_vec...
    %                     .*(gas_density_vec - melt_density)...
    %                     + gas_massfrac_vec*melt_density.*dgasdensity_dP_vec)...
    %                     - (gas_massfrac_vec - 1).*gas_density_vec.^2*dmeltdensity_dP)...
    %                     ./((-gas_massfrac_vec.*gas_density_vec + gas_massfrac_vec*melt_density + gas_density_vec).^2);
    ddensity_dP_vec(i) =  -((dgasmassfrac_dP_vec(i)/gas_density_vec(i)...
        - dgasmassfrac_dP_vec(i)/melt_density_vec(i) - (gas_massfrac_vec(i)*dgasdensity_dP_vec(i))...
        /gas_density_vec(i)^2 - ((1 - gas_massfrac_vec(i))*dmeltdensity_dP)...
        /melt_density_vec(i)^2)/(gas_massfrac_vec(i)/gas_density_vec(i)...
        + (1 - gas_massfrac_vec(i))/melt_density_vec(i))^2);    
    
end

soundspeed_vec = sqrt(1./ddensity_dP_vec);
bulkmod_vec = density_vec./ddensity_dP_vec;

phimax = max(porosity_vec);


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

% Ca_num = melt_viscosity*strain_rate*bubble_radius/surface_tension; %capillary number
%   Ca, which is function of fluid velocity and thus
%   of "fluid wave amplitude", which is not known prior to completing the
%   inversion...
% Cd_num = melt_viscosity*bubble_radius/surface_tension*wave_speed; %
%   dynamic capillary number, see Llewellin and Manga 2005
%   Cd_num obtained by assuming average fluid velocity is of the form uz = Acos(z - ct), 
%   which then gives avg|strain_rate|/avg|strain_rate_rate| = c
% magma_viscosity_ratio_bubble = fsolve(@(magma_viscosity_ratio)
% magma_viscosity_res_fun(magma_viscosity_ratio, Ca_num, porosity, porosity_cr), 1); %full magma viscosity model from Llewellin and Manga 2005 or Gonnermann and Manga 2012
% magma_viscosity_ratio_bubble = (1 - porosity)^(5/3); % function from Manga 2004 for when Cd,Ca >> 1
magma_viscosity_ratio_bubble = (1 - porosity_vec).^(-1); % function from Manga 2004 for when Cd,Ca << 1, which seems to be the case for reasonable ranges of halemaumau vlp wave parameters
% magma_viscosity_ratio_crystal = (1 - crystalvolfrac/crystalvolfrac_cr)^(-2.5); %Gonnermann and Manga 2012, Costa 2005 
mu_vec = magma_viscosity_ratio_bubble.*mu_melt_vec; 
% mu_vec = mu_melt_vec; %!!!ignoring bubbles!!!


%% convert to depth
Hlake_ind = find(zvec > Hlake,1); %index of closest depth to desired lake depth

restop_ind = find(zvec >= Hcolumn,1); %highest point within reservoir
if restop_ind == 1
    restop_ind = 2;
end
resbot_ind = find(zvec >= 2*Rres+Hcolumn,1) - 1; %lowest point within reservoir
P_restop = Pvec(restop_ind);

condtop_ind = find(zvec >= Hlake,1); %highest point within conduit
if condtop_ind == 1
    condtop_ind = 2;
end

%check for density stability
%                 density_unstable = any(mid_density_vec(2:end)<mid_density_vec(1:end-1)-stable_error_tol);
%check for thermodynamic stability (Kundu fluid mech ch 1.10)
%less strict than density stability, so redundant
potentialdensity = (movmean(density_vec(1:resbot_ind),11)...
    - movmean(density_vec(2:resbot_ind+1),11))./...
    (-dz)...
    + g*movmean(density_vec(1:resbot_ind),11).*...
    movmean(ddensity_dP_vec(1:resbot_ind),11);

thermo_unstable = any(potentialdensity < 0-stable_error_tol);


%% integrate over reservoir
delta_zvec_res = [zvec(restop_ind:resbot_ind);2*Rres+Hcolumn] - [Hcolumn;zvec(restop_ind:resbot_ind)];
mid_zvec_res = [Hcolumn;zvec(restop_ind:resbot_ind)] + delta_zvec_res/2; %mid points of each reservoir slice[cylinder]
%fractional location of midpoints between datapoints (accounts for reservoir bottom):
mid_frac_res = delta_zvec_res/2./(zvec(restop_ind:resbot_ind+1) - zvec(restop_ind-1:resbot_ind));
hor_radius_vec_res = real(sqrt(Rres^2 - (mid_zvec_res-Hcolumn-Rres).^2)); %radius at each midpoint
vol_vec_res = pi*hor_radius_vec_res.^2.*delta_zvec_res;%volume of each slice[cylinder]
mid_density_vec_res = density_vec(restop_ind-1:resbot_ind) ...
    + (density_vec(restop_ind:resbot_ind+1)...
    - density_vec(restop_ind-1:resbot_ind)).*mid_frac_res;
mid_ddensity_dP_vec_res = ddensity_dP_vec(restop_ind-1:resbot_ind) ...
    + (ddensity_dP_vec(restop_ind:resbot_ind+1)...
    - ddensity_dP_vec(restop_ind-1:resbot_ind)).*mid_frac_res;
% if ~thermo_unstable %|| ~density_unstable 
    %volume weighted sum over vectors
    density_res = sum(vol_vec_res.*mid_density_vec_res)/sum(vol_vec_res);
%                     '!!!need to account for changing layer boundaries with pressure!!!'
    ddensity_dP_res = sum(vol_vec_res.*mid_ddensity_dP_vec_res)/sum(vol_vec_res);
%                     ddensity_dP_new = sum(vol_vec./mid_density_vec.*mid_ddensity_dP_vec)...
%                         .*density/sum(vol_vec);
% else
%     density_res = NaN;
%     ddensity_dP_res = NaN;
% end


%% integrate over conduit
zvec_cond = [Hlake;zvec(Hlake_ind:restop_ind-1);Hcolumn];
% delta_zvec_cond = [zvec(Hlake_ind:restop_ind-1);Hcolumn] - [Hlake;zvec(Hlake_ind:restop_ind-1)];
delta_zvec_cond = zvec_cond(2:end)-zvec_cond(1:end-1);
mid_zvec_cond = (zvec_cond(2:end)+zvec_cond(1:end-1))/2;
%fractional location of midpoints between datapoints (accounts for reservoir bottom):
mid_frac_cond = delta_zvec_cond/2./(zvec(Hlake_ind:restop_ind) - zvec(Hlake_ind-1:restop_ind-1));
mid_Rvec_cond = Rcondtop - (Rcondtop - Rcondbot)/(Hcolumn-Hlake)*(mid_zvec_cond-Hlake);
Acondvec = pi*mid_Rvec_cond.^2;
% mid_Rvec_cond = Rcondtop; %!!!!no taper, for Chao's model!!!!
vol_vec_cond = pi*mid_Rvec_cond.^2.*delta_zvec_cond;%volume of each slice[cylinder]
mid_density_vec_cond = density_vec(Hlake_ind-1:restop_ind-1) ...
    + (density_vec(Hlake_ind:restop_ind)...
    - density_vec(Hlake_ind-1:restop_ind-1)).*mid_frac_cond;
d_density_dz_vec_cond = -(density_vec(Hlake_ind:restop_ind)...
    - density_vec(Hlake_ind-1:restop_ind-1))./...
    (zvec(Hlake_ind:restop_ind)...
    - zvec(Hlake_ind-1:restop_ind-1));
mid_ddensity_dP_vec_cond = ddensity_dP_vec(Hlake_ind-1:restop_ind-1) ...
    + (ddensity_dP_vec(Hlake_ind:restop_ind)...
    - ddensity_dP_vec(Hlake_ind-1:restop_ind-1)).*mid_frac_cond;
mid_soundspeed_vec_cond = soundspeed_vec(Hlake_ind-1:restop_ind-1) ...
    + (soundspeed_vec(Hlake_ind:restop_ind)...
    - soundspeed_vec(Hlake_ind-1:restop_ind-1)).*mid_frac_cond;
mid_bulkmod_vec_cond = bulkmod_vec(Hlake_ind-1:restop_ind-1) ...
    + (bulkmod_vec(Hlake_ind:restop_ind)...
    - bulkmod_vec(Hlake_ind-1:restop_ind-1)).*mid_frac_cond;
mid_mu_vec_cond = mu_vec(Hlake_ind-1:restop_ind-1) ...
    + (mu_vec(Hlake_ind:restop_ind)...
    - mu_vec(Hlake_ind-1:restop_ind-1)).*mid_frac_cond;
% if ~thermo_unstable %|| ~density_unstable 
    %volume weighted sum over vectors
    density_cond = sum(vol_vec_cond.*mid_density_vec_cond)/sum(vol_vec_cond);
%                     '!!!need to account for changing layer boundaries with pressure!!!'
    ddensity_dP_cond = sum(vol_vec_cond.*mid_ddensity_dP_vec_cond)/sum(vol_vec_cond);
%                     ddensity_dP_new = sum(vol_vec./mid_density_vec.*mid_ddensity_dP_vec)...
%                         .*density/sum(vol_vec);
    density_cond_top = mid_density_vec_cond(1);
    density_cond_bot = mid_density_vec_cond(end);
    mu_cond = sum(vol_vec_cond.*mid_mu_vec_cond)/sum(vol_vec_cond);
    
% else
%     density_cond = NaN;
%     ddensity_dP_cond = NaN;
%     density_cond_top = NaN;
%     density_cond_bot = NaN;
%     mu_cond = NaN;
% end


%% integrate over conduit+lake
zvec_condlake = [zvec(1:restop_ind-1);Hcolumn];
delta_zvec_condlake = zvec_condlake(2:end) - zvec_condlake(1:end-1);
mid_zvec_condlake = (zvec_condlake(2:end) + zvec_condlake(1:end-1))/2;
% delta_zvec_condlake = [zvec(2:restop_ind-1);Hcolumn] - zvec(1:restop_ind-1);
%fractional location of midpoints between datapoints (accounts for reservoir bottom):
mid_frac_condlake = delta_zvec_condlake/2./(zvec(2:restop_ind) - zvec(1:restop_ind-1));
mid_Rvec_condlake = zeros(size(mid_zvec_condlake)) + Rcondtop;

%!!!!!!conical lake!!!!!!
% mid_Rvec_condlake(mid_zvec_condlake < Hlake) = Rlake - ...
%     (Rlake-Rcondtop)/Hlake*mid_zvec_condlake(mid_zvec_condlake < Hlake);

%!!!!!!cylinder lake!!!!!!!!!
mid_Rvec_condlake(mid_zvec_condlake < Hlake) = Rlake; 
mid_Rvec_condlake(mid_zvec_condlake > Hlake) = Rcondtop - ...
    (Rcondtop - Rcondbot)/(Hcolumn-Hlake)*(mid_zvec_condlake(mid_zvec_condlake > Hlake)-Hlake);
Acondlakevec = pi*mid_Rvec_condlake.^2;
conddipvec = 90+zeros(size(mid_zvec_condlake));
conddipvec(mid_zvec_condlake > Hlake) = conddip;
mid_density_vec_condlake = density_vec(1:restop_ind-1) ...
    + (density_vec(2:restop_ind)...
    - density_vec(1:restop_ind-1)).*mid_frac_condlake;
mid_ddensity_dP_vec_condlake = ddensity_dP_vec(1:restop_ind-1) ...
    + (ddensity_dP_vec(2:restop_ind)...
    - ddensity_dP_vec(1:restop_ind-1)).*mid_frac_condlake;
mid_soundspeed_vec_condlake = soundspeed_vec(1:restop_ind-1) ...
    + (soundspeed_vec(2:restop_ind)...
    - soundspeed_vec(1:restop_ind-1)).*mid_frac_condlake;
mid_bulkmod_vec_condlake = bulkmod_vec(1:restop_ind-1) ...
    + (bulkmod_vec(2:restop_ind)...
    - bulkmod_vec(1:restop_ind-1)).*mid_frac_condlake;
d_density_dz_vec_condlake = -(density_vec(2:restop_ind)...
    - density_vec(1:restop_ind-1))./...
    (zvec(2:restop_ind)...
    - zvec(1:restop_ind-1));
mid_mu_vec_condlake = mu_vec(1:restop_ind-1) ...
    + (mu_vec(2:restop_ind)...
    - mu_vec(1:restop_ind-1)).*mid_frac_condlake;
% if ~thermo_unstable %|| ~density_unstable 
    density_condlake_top = mid_density_vec_condlake(1);
    
% else
%     density_condlake_top = NaN;
%     int_drhodz_over_A = NaN;
%     int_rho_over_A = NaN;
%     int_besselterm = NaN;
% end

%% integrate over lake
zvec_lake = [zvec(1:condtop_ind-1);Hlake];
delta_zvec_lake = zvec_lake(2:end) - zvec_lake(1:end-1);
mid_zvec_lake = (zvec_lake(2:end) + zvec_lake(1:end-1))/2;
% delta_zvec_lake = [zvec(2:condtop_ind-1);Hcolumn] - zvec(1:condtop_ind-1);
%fractional location of midpoints between datapoints (accounts for reservoir bottom):
mid_frac_lake = delta_zvec_lake/2./(zvec(2:condtop_ind) - zvec(1:condtop_ind-1));
mid_Rvec_lake = zeros(size(mid_zvec_lake)) + Rcondtop;

%!!!!!!conical lake!!!!!!
% mid_Rvec_lake(mid_zvec_lake < Hlake) = Rlake - ...
%     (Rlake-Rcondtop)/Hlake*mid_zvec_lake(mid_zvec_lake < Hlake);

%!!!!!!cylinder lake!!!!!!!!!
mid_Rvec_lake(mid_zvec_lake < Hlake) = Rlake; 

mid_Rvec_lake(mid_zvec_lake > Hlake) = Rcondtop - ...
    (Rcondtop - Rcondbot)/(Hcolumn-Hlake)*(mid_zvec_lake(mid_zvec_lake > Hlake)-Hlake);
%Alakevec = pi*mid_Rvec_lake.^2;
vol_vec_lake = pi*mid_Rvec_lake.^2.*delta_zvec_lake;%volume of each slice[cylinder]
mid_density_vec_lake = density_vec(1:condtop_ind-1) ...
    + (density_vec(2:condtop_ind)...
    - density_vec(1:condtop_ind-1)).*mid_frac_lake;
% d_density_dz_vec_lake = -(density_vec(2:condtop_ind)...
%     - density_vec(1:condtop_ind-1))./...
%     (zvec(2:condtop_ind)...
%     - zvec(1:condtop_ind-1));
mid_mu_vec_lake = mu_vec(1:condtop_ind-1) ...
    + (mu_vec(2:condtop_ind)...
    - mu_vec(1:condtop_ind-1)).*mid_frac_lake;
density_lake = sum(vol_vec_lake.*mid_density_vec_lake)/sum(vol_vec_lake);
mu_lake = sum(vol_vec_lake.*mid_mu_vec_lake)/sum(vol_vec_lake);


%% plot
if plotprofiles
  figure()
  %debug plot
  yticklist = (0:100:max(zvec)) + 0.1;
  yticklabelcell = cell(size(yticklist));
  for kk = 1:length(yticklabelcell)
      Pind = find(zvec >= yticklist(kk),1);
      yticklabelcell(kk) = strcat(num2str(yticklist(kk),'%.0f'),' m\newline',{' '},...
          num2str(Pvec(Pind)/1E6,'%.1f'),' MPa');
  end
  %ylimits = [min(zvec),max(zvec)];
  ylimits = [0,zvec(restop_ind)];
%   tiledlayout(1,3,'TileSpacing','compact') %figure()
figure()

%   nexttile%subplot(1,5,1)
subplot('position',[0.15,0.15,0.25,0.6])
  hold on
  plot(CO2_vec+H2O_vec,zvec,'k','LineWidth',2.5)
  plot(CO2_gas_vec+H2O_gas_vec,zvec,'k:','LineWidth',3)
%   plot(CO2_vec-CO2_gas_vec+H2O_vec-H2O_gas_vec,zvec,'r','LineWidth',3)
  plot(H2O_vec,zvec,'b','LineWidth',1)
  plot(H2O_gas_vec,zvec,'b:','LineWidth',2.4)
%   plot(H2O_vec-H2O_gas_vec,zvec,'r:','LineWidth',1.5)
  plot(CO2_vec,zvec,'r','LineWidth',1)
  plot(CO2_gas_vec,zvec,'r:','LineWidth',2.4)
%   plot(CO2_vec-CO2_gas_vec,zvec,'r--','LineWidth',1.5)
  xlabel('mass fraction (wt%)')
  %xlim([0,max(max(H2O_vec)*1.05,1E-8)])
  xlim([0,2.1])
%   ylabel('P (MPa), z (km)')
  yticks(yticklist)
  yticklabels(yticklabelcell)
  ylim(ylimits)
  grid on, box on
  legend({'total volatiles','total gas',...
      'H2O total','H2O gas',...
      'CO2 total','CO2 gas'},...
      'Location','northwest','AutoUpdate','off')
  set(gca,'ydir','reverse')
  yline(zvec(Hlake_ind));%,'k','conduit top','FontSize',14);
  yline(zvec(restop_ind));%,'k','conduit bottom','FontSize',14);
  yline(zvec(resbot_ind));%,'k','reservoir bottom','FontSize',14); 
  set(gca,'FontSize',14)

%   nexttile%subplot(1,5,2)
%   hold on
%   plot(CO2_vec,zvec,'k','LineWidth',1.5)
%   plot(CO2_gas_vec,zvec,':k','LineWidth',3)
%   xlabel('wt%')
%   xlim([0,max(max(CO2_vec)*1.05,1E-8)])
%   grid on, box on
%   legend('CO2 gas','CO2 total','Location','northwest','AutoUpdate','off')
%   set(gca,'ydir','reverse')
%   yticks(yticklist)
%   yticklabels([])
%   ylim(ylimits)
%   yline(zvec(Hlake_ind),'k','conduit top','FontSize',14);
%   yline(zvec(restop_ind),'k','conduit bottom','FontSize',14);
%   yline(zvec(resbot_ind),'k','reservoir bottom','FontSize',14); 
%   set(gca,'FontSize',14)

%   nexttile%subplot(1,5,3)
%   semilogx(porosity_vec,zvec,'Color','k','LineWidth',3)
%   hold on
%   semilogx(gas_massfrac_vec,zvec,':k','LineWidth',3)
%   xlim([1E-8,1])
%   xlabel('ratio')
%   grid on, box on
%   legend('porosity','gas\newlinemass fraction','Location','northwest','AutoUpdate','off')
%   set(gca,'ydir','reverse')
%   minlim = nanmin([porosity_vec(porosity_vec>0);...
%       gas_massfrac_vec(gas_massfrac_vec>0)])/2;
%   if length(minlim)==0
%       minlim = 0;
%   end
%   xlim([max(minlim,1E-8), 1])
%   yticks(yticklist)
%   yticklabels([])
%   ylim(ylimits)
% %   xticks(10.^(-20:0))
%   yline(zvec(Hlake_ind),'k','conduit top','FontSize',14);
%   yline(zvec(restop_ind),'k','conduit bottom','FontSize',14);
%   yline(zvec(resbot_ind),'k','reservoir bottom','FontSize',14); 
%   set(gca,'FontSize',14)

%   subplot(1,5,5)
%   plot(ddensity_dP_vec,zvec,'k','LineWidth',2)
%   hold on
%   plot(dgasdensity_dP_vec,zvec,':r','LineWidth',2)
%   plot(dmeltdensity_dP + 0*Pvec,zvec,'--c','LineWidth',1.5)
%   %semilogx(ddensity_dP*[1;1],...
%   %[P(Pcolumn_ind+top_ind);P(Pcolumn_ind+bot_ind+1)]/1E6,'*--g','LineWidth',0.5)
%   xlabel('\partialrho/\partialP')
%   grid on
%   legend('Bulk','Gas','Melt','Location','southeast','AutoUpdate','off')%'Integrated','Location','southeast')
%   set(gca,'ydir','reverse')
% %   xlim([dmeltdensity_dP/2,nanmax([dgasdensity_dP_vec;ddensity_dP_vec])*2])
%   yticks(yticklist)
%   yticklabels([])
%   ylim(ylimits)
%   xticks(10.^(-20:20))
%   yline(zvec(Hlake_ind),'k','conduit top');
%   yline(zvec(restop_ind),'k','conduit bottom');
%   yline(zvec(resbot_ind),'k','reservoir bottom'); 
%   set(gca,'XScale','Log')
%   set(gca,'FontSize',12)

%   nexttile%subplot(1,5,4)
  subplot('position',[0.425,0.15,0.15,0.6])
  plot(density_vec,zvec,'k','LineWidth',2.5)
  hold on
  plot(gas_density_vec,zvec,'k:','LineWidth',3)
  plot(melt_density_vec,zvec,'k--','LineWidth',3)
  %plot(density*[1;1],...
  %[P(Pcolumn_ind+top_ind);P(Pcolumn_ind+bot_ind+1)]/1E6,'*--g','LineWidth',0.5)
  xlabel({'magma density','\rho (kg/m^3)'})
  xlim([0,max(melt_density_vec)+10])
  grid on, box on
  legend('bulk','gas','melt','Location','northwest','AutoUpdate','off')%'Integrated','Location','southeast')
  set(gca,'ydir','reverse')
  yticks(yticklist)
  yticklabels([])
  ylim(ylimits)
  yline(zvec(Hlake_ind),'k-')%,'conduit top','FontSize',14);
  yline(zvec(restop_ind),'k-')%,'conduit bottom','FontSize',14);
  yline(zvec(resbot_ind),'k-')%,'reservoir bottom','FontSize',14);  
  set(gca,'FontSize',14)
  xticks([1200,2400])
  xlim([0,2800])
  
%   nexttile%subplot(1,5,5)
subplot('position',[0.6,0.15,0.15,0.6])
  plot(mu_vec,zvec,'k','LineWidth',2.5)
  hold on
  plot(mu_melt_vec,zvec,'k--','LineWidth',3)
  %semilogx(ddensity_dP*[1;1],...
  %[P(Pcolumn_ind+top_ind);P(Pcolumn_ind+bot_ind+1)]/1E6,'*--g','LineWidth',0.5)
  xlabel({'magma viscosity','\mu (Pas)'})
  grid on, box on
  legend('bulk','melt','Location','northwest','AutoUpdate','off')%'Integrated','Location','southeast')
  set(gca,'ydir','reverse')
%   xlim([dmeltdensity_dP/2,nanmax([dgasdensity_dP_vec;ddensity_dP_vec])*2])
  yticks(yticklist)
  yticklabels([])
  ylim(ylimits)
  xticks(10.^(-20:20))
  yline(zvec(Hlake_ind),'k-')%,'conduit top','FontSize',14);
  yline(zvec(restop_ind),'k-')%,'conduit bottom','FontSize',14);
  yline(zvec(resbot_ind),'k-')%,'reservoir bottom','FontSize',14); 
  set(gca,'XScale','Log')
  set(gca,'FontSize',14)

subplot('position',[0.8,0.15,0.15,0.6])
  plot(soundspeed_vec,zvec,'k','LineWidth',2.5)
  hold on
  %semilogx(ddensity_dP*[1;1],...
  %[P(Pcolumn_ind+top_ind);P(Pcolumn_ind+bot_ind+1)]/1E6,'*--g','LineWidth',0.5)
  xlabel({'sound speed','(m/s)'})
  grid on, box on
  legend('bulk','Location','northwest','AutoUpdate','off')%'Integrated','Location','southeast')
  set(gca,'ydir','reverse')
%   xlim([dmeltdensity_dP/2,nanmax([dgasdensity_dP_vec;ddensity_dP_vec])*2])
  yticks(yticklist)
  yticklabels([])
  ylim(ylimits)
  %xticks(10.^(-20:20))
  yline(zvec(Hlake_ind),'k-')%,'conduit top','FontSize',14);
  yline(zvec(restop_ind),'k-')%,'conduit bottom','FontSize',14);
  yline(zvec(resbot_ind),'k-')%,'reservoir bottom','FontSize',14); 
  %set(gca,'XScale','Log')
  set(gca,'FontSize',14)

%                   pause
end

%% output
output_struct = struct();
output_struct.density_res = density_res;
output_struct.ddensity_dP_res = ddensity_dP_res;
% output_struct.density_cond = density_cond;
% output_struct.ddensity_dP_cond = ddensity_dP_cond;
% output_struct.density_cond_top = density_cond_top;
% output_struct.mu_cond = mu_cond;
% output_struct.density_condlake_top = density_condlake_top;
% output_struct.density_cond_bot = density_cond_bot;

% output_struct.density_lake = density_lake;
% output_struct.mu_lake = mu_lake;

output_struct.mid_mu_vec_condlake = mid_mu_vec_condlake;
output_struct.mid_Rvec_condlake = mid_Rvec_condlake;
output_struct.Acondlakevec = Acondlakevec;
output_struct.mid_density_vec_condlake = mid_density_vec_condlake;
output_struct.mid_ddensity_dP_vec_condlake = mid_ddensity_dP_vec_condlake;
output_struct.mid_soundspeed_vec_condlake = mid_soundspeed_vec_condlake;
output_struct.mid_bulkmod_vec_condlake = mid_bulkmod_vec_condlake;
output_struct.delta_zvec_condlake = delta_zvec_condlake;
output_struct.mid_zvec_condlake = mid_zvec_condlake;
% output_struct.d_density_dz_vec_condlake = d_density_dz_vec_condlake;

output_struct.thermo_unstable = thermo_unstable;
% output_struct.phimax = phimax;
output_struct.P_restop = P_restop;
end
