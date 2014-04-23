#!/usr/bin/env python
# _Display_fit.py test_type line_from line_to P/M

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re

from ase.units import GPa
from ase.utils.eos import EquationOfState

#def centralfit(x, y):
#    list = []
#    x = np.array(x)
#    y = np.array(y)
#    for a in np.arange(160, 180, 0.1):
#        for c in np.arange(-14.1, -14, 0.001):
#            y_fit = a*(x**2)+c
#            ss_resid = np.sum((y - y_fit)**2)
#            list.append([ss_resid, a, c])
#    the_one = list[np.array(list)[:,0].argmin()]
#    y_fit = the_one[1]*(x**2)+the_one[2]
#    y_avg = np.sum(y)/len(y)
#    ss_resid = np.sum((y - y_fit)**2)
#    ss_total = np.sum((y-y_avg)**2)
#    r_squared = 1 - ss_resid/ss_total
#    return (the_one[1], the_one[2], r_squared, y_fit)
#

def polyfit(x, y, degree):
    popts = np.polyfit(x, y, degree, full=True)
    p = np.poly1d(popts[0])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 1000)
    y_fit = p(x_fit)
    y_avg = np.sum(y)/len(y)
    ss_total = np.sum((y-y_avg)**2)
    r_squared = 1 - popts[1]/ss_total
    return (p, r_squared[0], x_fit, y_fit)

#def murnaghan_eqn(V, V0, B0, B0_prime, E0):
#    return (E0 + B0/B0_prime * V * (1 + (V0/V)**B0_prime/(B0_prime - 1)) - V0 * B0/(B0_prime -1))
    
def murnaghan_eqn(V, V0, B0, B0_prime, E0):
    return (E0 + 9*V0*B0/16. * (((V0/V)**(2/3.) - 1)**3 * B0_prime + ((V0/V)**(2/3.) - 1)**2 * (6 - 4*(V0/V)**(2/3.))))

def murnaghan_fit(x, y):
    coeffs, pcov = curve_fit(murnaghan_eqn, x, y, [np.average(x), 2.5, 4, np.average(y)])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 1000)
    y_fit = murnaghan_eqn(x_fit, *coeffs)
    y_fit_eqlen = murnaghan_eqn(x, *coeffs)
    ss_resid = np.sum((y_fit_eqlen - y)**2)
    y_avg = np.sum(y)/len(y)
    ss_total = np.sum((y-y_avg)**2)
    r_squared = 1 - ss_resid/ss_total
    return (coeffs, r_squared, x_fit, y_fit)

test_type = sys.argv[1]
f = open(test_type+'_output.txt','rU')
file = f.readlines()
f.close()
line_from = int(sys.argv[2]); line_to = int(sys.argv[3])

if test_type == 'entest':
    ENCUT = []; dE1 = []; dE2 = []; dDE = []; 
    for i in file[line_from:line_to]:
        ENCUT.append(float(i.split()[0]))
        dE1.append(float(i.split()[2]))
        dE2.append(float(i.split()[4]))
        dDE.append(float(i.split()[6]))
    plt.plot(ENCUT, dE1, 'x', ENCUT, dE2, '*', ENCUT, dDE, 'o')
    plt.plot([ENCUT[0], ENCUT[-1]], [0.001, 0.001], 'k:')
    plt.grid()
    plt.xlabel('ENCUT (eV)')
    plt.ylabel('dDE (eV)')

elif test_type == 'kptest':
    nKP = []; dE1 = []; dE2 = []; dDE = []; 
    for i in file[line_from:line_to]:
        nKP.append(float(i.split()[0]))
        dE1.append(float(i.split()[2]))
        dE2.append(float(i.split()[4]))
        dDE.append(float(i.split()[6]))
    plt.plot(nKP, dE1, 'x', nKP, dE2, '*', nKP, dDE, 'o')
    plt.plot([nKP[0], nKP[-1]], [0.001, 0.001], 'k:')
    plt.grid()
    plt.xlabel('nKP')
    plt.ylabel('dDE (eV)')

elif test_type == 'lctest':
    scaling_factor = []; volume = []; energy = []
    for i in file[line_from:line_to]:
        scaling_factor.append(float(i.split()[0]))
        volume.append(float(i.split()[1]))
        energy.append(float(i.split()[2]))
    scaling_factor = np.array(scaling_factor)
    volume = np.array(volume)
    energy = np.array(energy)
    V_a_conversion_multiplier = scaling_factor[0]**3/volume[0]

    plt.plot(volume, energy, 'o', label="Original data")

    # fitting the polynomial
#    (p_vol, r_squared, volume_fit, energy_fit) = polyfit(volume, energy, 3)
#    volume_eqlbrm_by_polynomial = volume_fit[energy_fit.argmin()]
#    scaling_factor_eqlbrm_by_polynomial = (volume_eqlbrm_by_polynomial*V_a_conversion_multiplier)**(1/3.)
    
    # fitting the Birch-Murnaghan equation of state
    (coeffs_vol_M, r_squared_M, volume_fit_M, energy_fit_M) = murnaghan_fit(volume, energy)
    scaling_factor_eqlbrm = (coeffs_vol_M[0]*V_a_conversion_multiplier)**(1/3.)

    # fitting B-M eos with ase module
#    eos = EquationOfState(volume, energy)
#    volume_eqlbrm_by_ase, energy_eqlbrm_by_ase, B_by_ase = eos.fit()

    # plotting the Birch-Murnaghan equation of state
    plt.plot(volume_fit_M, energy_fit_M, '-', label="Birch-Murnaghan eqn of state")
    result_str = "R-squared is %f" % r_squared_M
    plt.text(volume_fit_M[len(volume_fit_M)/4], energy_fit_M[60], result_str)

    # standrad output, directed to files by the bash script calling this python script
    print "%s" % result_str
#    print("Equilibrium scaling factor is {0} , or {1} (polynomial)".format(scaling_factor_eqlbrm, scaling_factor_eqlbrm_by_polynomial))
    print("Equilibrium scaling factor is {0}".format(scaling_factor_eqlbrm))
    if scaling_factor_eqlbrm <= scaling_factor[0] or scaling_factor_eqlbrm >= scaling_factor[-1]:
        print("!Equilibrium point is out of the considered range!")
    else:
        print("V0 = %f\nB0 = %f\nB0' = %f" % (coeffs_vol_M[0], coeffs_vol_M[1] * 160.2, coeffs_vol_M[2]))
#        print("Total energy is {0} , or {1} (minimum of polynomial)".format(coeffs_vol_M[3], energy_fit.min()))
        print("Total energy is {0}".format(coeffs_vol_M[3]))
#        print "\nEquation of state parameters by ase"
#        print "----------------------------"
#        print "V0: %8.4f A^3" % (volume_eqlbrm_by_ase)
#        print "B: %8.2f GPa" % (B_by_ase/GPa)
#        print "E0: %8.6f eV" % (energy_eqlbrm_by_ase)

        np.savetxt('eosfit_data.dat', np.column_stack((volume_fit_M, energy_fit_M)), '%.6f', '\t')

    plt.xlabel(r'Volume ($\AA^{3}$)')
    plt.ylabel('E (eV)')
    plt.legend(loc=0)
    np.savetxt('orig_data.dat', np.column_stack((volume, energy)), '%.6f', '\t')

elif test_type == 'rttest':
    ratio = []; volume = []; energy = []
    for i in file[line_from:line_to]:
        ratio.append(float(i.split()[0]))
        volume.append(float(i.split()[1]))
        energy.append(float(i.split()[2]))
    ratio = np.array(ratio)
    volume = np.array(volume)
    energy = np.array(energy)
    plt.plot(ratio, energy, 'o')
    (p_ratio, r_squared, ratio_fit, energy_fit) = polyfit(ratio, energy, 4)
    ratio_eqlbrm = ratio_fit[energy_fit.argmin()]
    plt.plot(ratio_fit, energy_fit, '-')

    print("R-squared is {0}\nEquilibrium ratio is {1}\nThe polynomial is\n{2}".format(r_squared, ratio_eqlbrm, p_ratio))
    print("Minimal total energy is %f" % energy_fit.min())
    np.savetxt('polyfit_data.dat', np.column_stack((ratio_fit, energy_fit)), '%.6f', '\t')
    plt.xlabel('raito')
    plt.ylabel('E (eV)')
    np.savetxt('orig_data.dat', np.column_stack((ratio, energy)), '%.6f', '\t')

elif test_type == 'agltest':
    angle = []; energy = []
    for i in file[line_from:line_to]:
        angle.append(float(i.split()[0]))
        energy.append(float(i.split()[1]))
    angle = np.array(angle)
    energy = np.array(energy)
    plt.plot(angle, energy, 'o')
    plt.xlabel('angle')
    plt.ylabel('E (eV)')
    np.savetxt(test_type+'_orig_data.dat', np.column_stack((angle, energy)), '%.6f', '\t')

elif re.search('.*c[1-9][1-9].*', test_type) or re.search('A.*', test_type):     # meaning elastic const.
    delta = []; energy = []
    for i in file[line_from:line_to]:
        delta.append(float(i.split()[0]))
        energy.append(float(i.split()[1]))
        
    if re.search('.*c[1-9][1-9].*', test_type):
        (p, r_squared, delta_fit, energy_fit) = polyfit(delta, energy, 2)
        plt.plot(delta, energy, 'o', delta_fit, energy_fit, '-')
        result_str = "E = %f x^2 + (%f) x + (%f)\nR-squared is %f" % (p[2], p[1], p[0], r_squared)
    elif re.search('A.*', test_type):
        (p, r_squared, delta_fit, energy_fit) = polyfit(delta, energy, 3)
        plt.plot(delta, energy, 'o', delta_fit, energy_fit, '-')
        result_str = "E = %f x^3 + (%f) x^2 + (%f) x + (%f)\nR-squared is %f" % (p[3], p[2], p[1], p[0], r_squared)
        
    plt.text(delta_fit[len(delta_fit)/4], energy_fit[6], result_str)
    plt.xlabel('Delta (ratio)')
    plt.ylabel('E (eV)')
    print "Fitting result:", result_str
    

plt.savefig('fit_curve.png')
#plt.show()
