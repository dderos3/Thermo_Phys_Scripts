import numpy as np
import matplotlib as plt
"""
TODO: 
        - Add functions to plot the distribution
        - Add error handling for n < 0
        - Add separate file for constants and useful tables
        - Add dictonary for elements and their atomic masses
        - Test worst cases to see if float64 casting is really nessesary to prevent overflows
"""
# Constants
K_J = np.float64(1.38*pow(10, -23))         # Boltzman constant 'K' in J/k
K_EV = np.float64(8.617*pow(10, -5))        # Boltzman constant 'K' in ev/k
AMU_KG = np.float64(1.66*pow(10, -27))      # Atomic mass unit Kg conversion scalar 

# System parameters
n = 5                                       # number of energy levels
Temp = 300                                  # Temperature in Kelvin
epsilon = np.float64(pow(10, -20))          # Energy scale
atm_mass = 32                               # Atomic mass of the element


# Returns the probability of the power level
def CalcEnergyLvlProb(n, T, epsilon, power_level):
    sum = 0
    for i in range(0, n):
        sum = sum + np.float64(pow(np.e, ((-1*i*epsilon) / (K_J * T))))
    p_power_level = np.float64(pow(np.e, ((-1*power_level*epsilon) / (K_J * T))))
    return np.float64((p_power_level / sum))

# Returns the average energy in the system
def CalcAvgEnergy():
    sum = 0
    for i in range(0, n):
        sum = sum + (CalcEnergyLvlProb(n, Temp, epsilon, i) * i)
    return sum

# Returns the max velocity of a mol. at a given temp
def CalcVMax(atm_mass, T):
    mass_kg = atm_mass * AMU_KG
    return (np.sqrt((2*K_J*T) / mass_kg))

# Returns the average velocity of a mol. at a given temp
def CalcVAvg(atm_mass, T):
    mass_kg = atm_mass * AMU_KG
    return (np.sqrt((8*K_J*T) / (np.pi * mass_kg)))


