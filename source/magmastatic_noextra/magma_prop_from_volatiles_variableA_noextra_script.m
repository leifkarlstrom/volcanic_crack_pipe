%script to calculate magmastatic properties

conddip = 90;
Rres = 1250;
Zrescentroid = 1100-1940;
Zlaketop = 900;
Zlakebot = 700;
Zcondbot = Zrescentroid + Rres;
Hcolumn = Zlaketop - Zcondbot;
Hlake = Zlaketop - Zlakebot;
Rcondtop = 15;
Rcondbot = 15;
Rlake = 100;

TK_condbot = 1150 + 273.15;
TK_condtop = 1150 + 273.15;
TK_lakebot = 1150 + 273.15;
TK_laketop = 1150 + 273.15;

H2O_laketop_wtp = 1.5;
H2O_lakebot_wtp = 1.5;
H2O_condtop_wtp = 1.5;
H2O_condbot_wtp = 0.5;

wtp_CO2_to_ppm = 1/(7.0778e-05);
CO2_laketop_ppm = 1/3*H2O_laketop_wtp*wtp_CO2_to_ppm;
CO2_lakebot_ppm = 1/3*H2O_lakebot_wtp*wtp_CO2_to_ppm;
CO2_condtop_ppm = 1/3*H2O_condtop_wtp*wtp_CO2_to_ppm;
CO2_condbot_ppm = 1/3*H2O_condbot_wtp*wtp_CO2_to_ppm;

dz = 0.4; %depth increment for integration, recommend 0.4 or less

plotprofiles = true;

%http://calcul-isto.cnrs-orleans.fr/thermodynamics/applications/#/script/js_compute~H2O-CO2%20systems~H2O%20and%20CO2%20solubility.js
CO2solppm_table = readmatrix('new_Values_CO2solubility_ppm_comb.csv');
H2Osolwtp_table = readmatrix('new_Values_H2Osolubility_wtp_comb.csv');


magma_prop_struct = magma_prop_from_volatiles_variableA_noextra(conddip,...
    H2O_condtop_wtp,H2O_condbot_wtp,CO2_condtop_ppm,CO2_condbot_ppm,...
    H2O_laketop_wtp,H2O_lakebot_wtp,CO2_laketop_ppm,CO2_lakebot_ppm,...
    Hlake,Hcolumn,Rres,Rcondtop,Rcondbot,Rlake,...
    TK_condbot,TK_condtop,TK_lakebot,TK_laketop,...
    plotprofiles, CO2solppm_table, H2Osolwtp_table, dz);