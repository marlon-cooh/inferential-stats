import numpy as np

def flow_rate(u, c, area):
    # Returns flow rate based on initial reagent concentration, velocity flow, and reactor area
    # Area is the parameter to iterate for this analysis.
    return (u * c * area)

def mass_balance(area_r, r_ammonia, eta, u, c_0):
    # Define flow rate
    f_n = flow_rate(u, c_0, area_r)
    
    # Define ordinary differential equation
    return (r_ammonia * eta * area_r / 2 * f_n)

def energy_balance(area_r, r_ammonia, enthalpy, eta, u, c_0, cp_mix):
    # Define flow rate
    f_t = flow_rate(u, c_0, area_r)
    
    # Define ordinary differential equation
    return (-enthalpy * eta * r_ammonia / f_t * cp_mix)

def ammonia_rate(k, k_eq, a_n2, a_h2, a_nh3, alpha):
    # Returns ammonia rate
    r_ammonia = 2 * k * (k_eq**2 * a_n2 * (a_h2**3 / a_nh3**2)**alpha - (a_nh3**2 / a_h2**3)**(1-alpha))

    return r_ammonia

def molar_fraction(i, reagents):
    # Returns molar fraction
    y = reagents[i] / sum(reagents.values())
    return y

def fugacities(T, P): # Units probably do not match these empirical formulas.
    # Returns fugacity of hydrogen
    phi_h2 = np.exp(np.exp(-3.7301 * T**0.125 + 0.543) * P - np.exp(-0.1374 * T**0.5 - 16.99) * P**2 + 314 * np.exp(-0.02291 * T - 5.943) * np.exp(-P/300))
    
    # Returns fugacity of nitrogen
    phi_n2 = 0.94231827 + 0.2129547 * 10**(-3) * T + 28678509 * 10**(-3) * P - 0.280826 * 10**(-6) * T**2 + 0.4886308 * 10**(-6) * P**2

    # Returns fugacity of ammonia
    phi_nh3 = 0.1327886 + 0.3038548 * 10**(-2) * T - 0.4397572 * 10**(-3) * P - 0.1233954 * 10**(-5) * T**2 + 0.2851314 * 10**(-6) * P**2
    
    return (phi_h2, phi_n2, phi_nh3)

def activities(T, P, reagents):
    # Returns activities based on fugacities, operational pressure and molar fraction
    a_h2 = fugacities(T, P)[0] * molar_fraction('H_2', reagents)
    a_n2 = fugacities(T, P)[1] * molar_fraction('N_2', reagents)
    a_nh3 = fugacities(T, P)[2] * molar_fraction('NH_3', reagents)
    return (a_h2, a_n2, a_nh3)

def arrhenius(T):
    # Returns velocity rate constant according to Arrhenius
    k_0 = 8.85e14
    e_a = 170550 # kJ/kmol
    R = 8.314 # kJ/(kmol*K)
    k = k_0 * np.exp(-e_a/R*T)
    return k

def equilibrium_constant(T_exp):
    # Returns equilibria constant
    # k_eq = -2.792312 * np.log(T) - 5.527463 * 10**(-5) * T + 1.837742 * 10**(-5) * T**2 + 2002.4 / T + 2.799 # No convergence
    
    # Standard values
    delta_H_std = -92.22 # kJ/mol
    delta_gibbs_std = -32.8 # kJ/mol*K
    R = 8.314 # J/(mol*K)
    T = 298.15 # K
    k_eq_std = np.exp((-delta_gibbs_std*1000) / (T*R)) # Equilibria constant at SATP

    # Equilibria constant based on Van't Hoff equation
    v_hoff = (-1000*delta_H_std / R) * ((1/T) - (1/T_exp))
    k_eq = k_eq_std * np.exp(v_hoff)
    return k_eq

def catalyst_eff(T, X=1):

    coefficients = {
    "b0": -8.21,
    "b1": 0.038,
    "b2": 6.19,
    "b3": -5.35,
    "b4": -20.8,
    "b5": 2.38e-8,  # Using scientific notation for 2.38 × 10^-8
    "b6": 27.9
    }

    eta = (
    coefficients["b0"] +
    coefficients["b1"] * T +
    coefficients["b2"] * X +
    coefficients["b3"] * T**2 +
    coefficients["b4"] * X**2 +
    coefficients["b5"] * T**3 +
    coefficients["b6"] * X**3
    )

    return eta

def run():
    pass

if __name__ == "__main__":
    run()