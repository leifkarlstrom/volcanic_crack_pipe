import pandas as pd
import PetThermoTools as ptt
import Thermobar as pt
import numpy as np
import sys
import warnings


## This is a PYTHON module ##
## If an error occurs calling this module, first check if an error is occuring around line 69 with type(P_Pa)

sys.path.append(r"/Users/kwilde/Documents/phd_research/alphamelts/alphamelts_py/alphamelts-py-2.3.1-macos-arm64")

def run_PTT(n_tot, H2Ofrac, P_Pa, T_K):
    """
    Compute gas and melt properties based on given parameters.

    Args:
        n_tot (float): Total amount of H2O and CO2 (mass fraction).H2Ofrac (float): Fraction of n_tot that is H2O.
        P_Pa (float, int, or array): Pressure in Pa.
        T_K (float): Temperature in Kelvin.

    Returns:
        dict: Results including gas density, melt density, CO2 percent, H2O percent, and bulk viscosity.
    """
    warnings.filterwarnings("ignore", category=FutureWarning, message=".*chained assignment.*")

    # Define initial composition (non-H2O and non-CO2 components)
    # anhydrous composition
    # adds to 100 %
    wtp_sum_nonH2O = 50.6 + 2.3 + 12.9 + 10.56 + 0.87 + 0.19 + 8.95 + 10.6 + 2.19 + 0.38 + 0.23 + 0
    wtp_scale = 100 / wtp_sum_nonH2O  # Adjust so that % adds to 100
    
    wtp_comps = {
        'SiO2': wtp_scale * 50.6,
        'TiO2': wtp_scale * 2.3,
        'Al2O3': wtp_scale * 12.9,
        'FeO': wtp_scale * 10.56,
        'Fe2O3': wtp_scale * 0.87,
        'MnO': wtp_scale * 0.19,
        'MgO': wtp_scale * 8.95,
        'CaO': wtp_scale * 10.6,
        'Na2O': wtp_scale * 2.19,
        'K2O': wtp_scale * 0.38,
        'P2O5': wtp_scale * 0.23,
        'F2O_1': wtp_scale * 0  # Assuming
    }



    # calculate wt% of H2O and CO2
    wtp_volatiles_tot = n_tot * 100
    wtp_H2O = round(wtp_volatiles_tot * H2Ofrac, 5)
    wtp_CO2 = round(wtp_volatiles_tot - wtp_H2O, 5)
    
    # scale wt% H2O and CO2 
    wtp_H2O_scaled = round((wtp_H2O * (100 + wtp_volatiles_tot)) / 100, 5)
    wtp_CO2_scaled = round((wtp_CO2 * (100 + wtp_volatiles_tot)) / 100, 5)

    
    # Add H2O and CO2 to the composition
    wtp_comps['H2O'] = wtp_H2O_scaled
    wtp_comps['CO2'] = wtp_CO2_scaled

    

    # needs to be input in deg C
    T_C = T_K - 273.15


    # this is completely necessary to avoid a slew of issues
    # if you are getting cryptic errors, start here...
    if isinstance(P_Pa, float) or isinstance(P_Pa, int):
        # Convert pressures 
        P_bar = P_Pa / 1E5
        
        # Create a single-row DataFrame
        input_df = pd.DataFrame([{**{'P_bar': P_bar, 'T_C': T_C}, **wtp_comps}])
        
        # Run equilibrium calculation
        out = ptt.equilibrate_multi(Model = "MELTSv1.2.0", bulk = wtp_comps, P_bar = input_df['P_bar'], 
                                    T_C = input_df['T_C'], fO2_buffer = "FMQ", Suppress = "All")
    
    else:
        
        # Convert pressures
        P_bar = [p / 1E5 for p in P_Pa]
        
        # Create dataframe
        input_df = pd.DataFrame([{**{'P_bar': p, 'T_C': T_C}, **wtp_comps} for p in P_bar])
        
         # Run equilibrium calculation
        out = ptt.equilibrate_multi(Model = "MELTSv1.2.0", bulk = input_df, P_bar = input_df['P_bar'], T_C = input_df['T_C'], fO2_buffer = "FMQ", Suppress = "All")

    
    # Calculate viscosity
    Vis = pt.calculate_viscosity_giordano_2008(liq_comps=out, T=out['T_C'] + 273.15)

    # Replace NaN values in mass_fluid1 with 0
    out["mass_fluid1"] = out["mass_fluid1"].fillna(0)
    out["v_fluid1"] = out["v_fluid1"].fillna(0)  

    # Replace NaN values with zeroes for porosity calculation
    out["rho_fluid1"] = out["rho_fluid1"].fillna(0) 
    
    # Calculate bulk density
    out["rho_bulk"] = (out["mass_Liq"] + out["mass_fluid1"]) / (out["v_Liq"] + out["v_fluid1"])
    
    # Calculate ngas (gas mass fraction)
    out["n_gas"] = out["mass_fluid1"] / (out["mass_Liq"] + out["mass_fluid1"])
    
    # Calculate volume percent gas
    out["vpercent_gas"] = out["v_fluid1"] / (out["v_Liq"] + out["v_fluid1"])

     # Extract required outputs
     # Convert to kg/m3 from g/cm3
    results = {
        'gas_density': out['rho_fluid1'].values * 1000,
        'melt_density': out['rho_Liq'].values * 1000,
        'melt_viscosity': Vis['n_melt'].values,
        'bulk_density': out["rho_bulk"].values * 1000,
        'n_gas':  out["n_gas"].values
    }

    # Need this to make it Matlab compatible...
    results = {key: value.tolist() for key, value in results.items()}  

    return results