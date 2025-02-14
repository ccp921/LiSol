'''
Python script as the interface between Cantera Python interface and Fortran program LiSol
Input: the mechanism filename in fort2py.dat
Output: 1 - Number of each atom type in the corresponding species --> *_spec.dat
        2 - Thermo data of each species, including Tlow, Thigh, Tmid, and 7 coeffs in highT range and lowT range of NASA7 (only for NASA7 now) --> *_thermo.dat
        3 - Trans data of each species, including geometry(0-single atom, 1-linear, 2-nonlinear), LJ well depth(K), LJ diameter(A), dipole, polarizability, and rotational relaxation --> *_trans.dat
Author: Chongpeng Chen
Email: chongpengchen@buaa.edu.cn
Date: 2025/01
'''
import sys
sys.path.insert(0, '/home/ccp921/cantera-opt/install/lib/python3.8/site-packages')  # specify your Cantera location here, used for different versions
import cantera as ct
import numpy as np
import warnings

# constans
boltzmann = ct.boltzmann  # Boltzmann constant, J/K
m2A = 1e10  # meter --> Angstrom
LightSpeed = ct.light_speed  # light speed, m/s
cal2J = 4.184
# reaction types that are supported now
ELEMENTARY_RXN = 100
THREE_BODY_RXN = 200
FALLOFF_RXN = 400
SIMPLE_FALLOFF = 410
TROE_FALLOFF = 420
SRI_FALLOFF_3 = 433
SRI_FALLOFF_5 = 435
PLOG_RXN = 500
non = 0

# same as LiSol
file_in = './fort2py.dat'
inter_path = './exchange_data'
spec_file_suffix = '_spec.dat'
thermo_file_suffix = '_thermo.dat'
trans_file_suffix = '_trans.dat'
Total_reac_type = 'reac_type.dat'
reac_file_suffix = '_reac.dat'

# get the mechanism file name with absolute path
fi = open(file_in, 'r')
lines = fi.readlines()
mech = lines[0][0:-1]
# print(mech)
fi.close()

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    gas = ct.Solution(mech)

#****************************************
# species properties
elem_name = gas.element_names
nelem = gas.n_elements
for species in gas.species():
    # element atom number
    composition = species.composition
    elem_number = np.zeros((nelem,))
    spec_file = inter_path+'/'+str(species.name)+spec_file_suffix
    fspec = open(spec_file, 'w')
    for i in range(nelem):
        for element, count in composition.items():
            if str(element) == elem_name[i]:
                elem_number[i] = count
                break
        fspec.write('{:.0f}\n'.format(elem_number[i])) # in the order of elements in the mechanism
    fspec.close()

    # species thermodynamic data
    thermo = species.thermo
    thermo_file = inter_path+'/'+str(species.name)+thermo_file_suffix
    fthermo = open(thermo_file, 'w')
    fthermo.write('{:.2f}\n'.format(thermo.min_temp))      # low T
    fthermo.write('{:.2f}\n'.format(thermo.coeffs[0]))     # middle T
    fthermo.write('{:.2f}\n'.format(thermo.max_temp))      # high T
    for i in range(8,15):
        fthermo.write('{:.16e} '.format(thermo.coeffs[i]))  # coeffs in low T range
    fthermo.write('\n')
    for i in range(1,8):
        fthermo.write('{:.16e} '.format(thermo.coeffs[i]))  # coeffs in high T range
    fthermo.close()

    # species transport data
    trans = species.transport
    trans_file = inter_path+'/'+str(species.name)+trans_file_suffix
    ftrans = open(trans_file, 'w')
    if(str(trans.geometry) == 'atom'):
        ftrans.write('{:.0f}\n'.format(0)) # single atom
    elif(str(trans.geometry) == 'linear'):
        ftrans.write('{:.0f}\n'.format(1)) # linear molecule
    elif(str(trans.geometry) == 'nonlinear'):
        ftrans.write('{:.0f}\n'.format(2)) # nonlinear molecule
    ftrans.write('{:.3f}\n'.format(trans.well_depth / boltzmann))  # LJ well depth, J --> K
    ftrans.write('{:.3f}\n'.format(trans.diameter * m2A))  # LJ diameter, meter --> Angstrom
    ftrans.write('{:.3f}\n'.format(trans.dipole * LightSpeed / 1e-21))  # dipole moment, Coulomb*m --> Debye
    ftrans.write('{:.3f}\n'.format(trans.polarizability * (m2A**3)))  # polarizability, m^3 --> Angstrom^3
    ftrans.write('{:.3f}\n'.format(trans.rotational_relaxation))  # rotational_relaxation
    ftrans.close()
#****************************************
#****************************************
# reaction properties
max_Pnum = 0
reac_total_file = inter_path+'/'+Total_reac_type
freac_t = open(reac_total_file, 'w')
for i, reaction in enumerate(gas.reactions()):
    reac_file = inter_path+'/'+str(i+1)+reac_file_suffix
    freac = open(reac_file, 'w')
    # 2. three-body reaction
    if isinstance(reaction, ct.ThreeBodyReaction):
        if reaction.efficiencies:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(THREE_BODY_RXN, len(reaction.efficiencies.items()), non))
        else:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(THREE_BODY_RXN, non, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.rate.pre_exponential_factor,\
                                                    reaction.rate.temperature_exponent,\
                                                    reaction.rate.activation_energy)) # J/kmol  
        if reaction.efficiencies:
            for species, efficiency in reaction.efficiencies.items():
                freac.write(f"{species}\n")
                freac.write(f"{efficiency}\n")
    # 3. simple falloff reaction
    elif isinstance(reaction, ct.FalloffReaction) and ('Simple' == reaction.falloff.type):
        if reaction.efficiencies:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SIMPLE_FALLOFF, len(reaction.efficiencies.items())))
        else:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SIMPLE_FALLOFF, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.low_rate.pre_exponential_factor,\
                                                    reaction.low_rate.temperature_exponent,\
                                                    reaction.low_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.high_rate.pre_exponential_factor,\
                                                    reaction.high_rate.temperature_exponent,\
                                                    reaction.high_rate.activation_energy))
        if reaction.efficiencies:
            for species, efficiency in reaction.efficiencies.items():
                freac.write(f"{species}\n")
                freac.write(f"{efficiency}\n")
    # 4. TROE falloff reaction
    elif isinstance(reaction, ct.FalloffReaction) and ('Troe' == reaction.falloff.type):
        if reaction.efficiencies:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, TROE_FALLOFF, len(reaction.efficiencies.items())))
        else:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, TROE_FALLOFF, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.low_rate.pre_exponential_factor,\
                                                    reaction.low_rate.temperature_exponent,\
                                                    reaction.low_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.high_rate.pre_exponential_factor,\
                                                    reaction.high_rate.temperature_exponent,\
                                                    reaction.high_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e} {:.16e}\n'.format(*np.array(reaction.falloff.parameters)))
        if reaction.efficiencies:
            for species, efficiency in reaction.efficiencies.items():
                freac.write(f"{species}\n")
                freac.write(f"{efficiency}\n")
    # 5. SRI-3 falloff reaction
    elif isinstance(reaction, ct.FalloffReaction) and (len(np.array(reaction.falloff.parameters))==3) and ('SRI' == reaction.falloff.type):
        if reaction.efficiencies:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SRI_FALLOFF_3, len(reaction.efficiencies.items())))
        else:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SRI_FALLOFF_3, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.low_rate.pre_exponential_factor,\
                                                    reaction.low_rate.temperature_exponent,\
                                                    reaction.low_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.high_rate.pre_exponential_factor,\
                                                    reaction.high_rate.temperature_exponent,\
                                                    reaction.high_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(*np.array(reaction.falloff.parameters)))
        if reaction.efficiencies:
            for species, efficiency in reaction.efficiencies.items():
                freac.write(f"{species}\n")
                freac.write(f"{efficiency}\n")
    # 6. SRI-5 falloff reaction
    elif isinstance(reaction, ct.FalloffReaction) and (len(np.array(reaction.falloff.parameters))==5) and ('SRI' == reaction.falloff.type):
        if reaction.efficiencies:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SRI_FALLOFF_5, len(reaction.efficiencies.items())))
        else:
            freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(FALLOFF_RXN, SRI_FALLOFF_5, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.low_rate.pre_exponential_factor,\
                                                    reaction.low_rate.temperature_exponent,\
                                                    reaction.low_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.high_rate.pre_exponential_factor,\
                                                    reaction.high_rate.temperature_exponent,\
                                                    reaction.high_rate.activation_energy))
        freac.write('{:.16e} {:.16e} {:.16e} {:.16e} {:.16e}\n'.format(*np.array(reaction.falloff.parameters)))
        if reaction.efficiencies:
            for species, efficiency in reaction.efficiencies.items():
                freac.write(f"{species}\n")
                freac.write(f"{efficiency}\n")
    # 7. PLOG reactions
    elif isinstance(reaction, ct.PlogReaction):
        if len(reaction.rates) > max_Pnum:
            max_Pnum = len(reaction.rates)
        freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(PLOG_RXN, len(reaction.rates), non))
        for rate in reaction.rates:
            freac.write('{:.16e} {:.16e} {:.16e} {:.16e}\n'.format(rate[0], \
                                                               rate[1].pre_exponential_factor, \
                                                               rate[1].temperature_exponent, \
                                                               rate[1].activation_energy))
    # 1. elementary reaction with Arrhenius rate
    elif isinstance(reaction, ct.ElementaryReaction):
        freac_t.write('{:.0f} {:.0f} {:.0f}\n'.format(ELEMENTARY_RXN, non, non))
        freac.write('{:.16e} {:.16e} {:.16e}\n'.format(reaction.rate.pre_exponential_factor,\
                                                    reaction.rate.temperature_exponent,\
                                                    reaction.rate.activation_energy))
    else:
        print('Unsupported reaction type for reaction {:.0f}!'.format(i+1))
        sys.exit(0)
    freac.close()
freac_t.write('{:.0f}\n'.format(max_Pnum))
freac_t.close()
#****************************************