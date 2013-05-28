#!/usr/local/python/2.7.1/bin/python
# _Display_fit.py TestType LineFrom LineTo

import sys
import numpy as np
import matplotlib.pyplot as plt
import re

def centralfit(x, y):
    list = []
    x = np.array(x)
    y = np.array(y)
    for a in np.arange(160, 180, 0.1):
        for c in np.arange(-14.1, -14, 0.001):
            y_fit = a*(x**2)+c
            ss_resid = np.sum((y - y_fit)**2)
            list.append([ss_resid, a, c])
    the_one = list[np.array(list)[:,0].argmin()]
    y_fit = the_one[1]*(x**2)+the_one[2]
    y_avg = np.sum(y)/len(y)
    ss_resid = np.sum((y - y_fit)**2)
    ss_total = np.sum((y-y_avg)**2)
    r_squared = 1 - ss_resid/ss_total
    return (the_one[1], the_one[2], r_squared, y_fit)

def polyfit(x, y, degree):
    coefficients = np.polyfit(x, y, degree, full=True)
    p = np.poly1d(coefficients[0])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 100)
    y_fit = p(x_fit)
    y_avg = np.sum(y)/len(y)
    ss_total = np.sum((y-y_avg)**2)
    r_squared = 1 - coefficients[1]/ss_total
    return (p, r_squared[0], x_fit, y_fit)

test_type = sys.argv[1]
f = open(test_type+'_output.txt','rU')
file = f.readlines()
f.close()
line_from = int(sys.argv[2]); line_to = int(sys.argv[3])

if test_type == 'entest':
    ENCUT = []; energy = []
    for i in file[line_from:line_to]:
        ENCUT.append(float(i.split()[0]))
        energy.append(float(i.split()[3]))
    plt.plot(ENCUT, energy, '.')
    plt.xlabel('ENCUT (eV)')
    plt.ylabel('E (eV)')

elif test_type == 'kptest':
    nKP = []; energy = []
    for i in file[line_from:line_to]:
        nKP.append(float(i.split()[0]))
        energy.append(float(i.split()[3]))
    plt.plot(nKP, energy, '.')
    plt.xlabel('nKP')
    plt.ylabel('E (eV)')

elif test_type == 'lctest':
    scaling_factor = []; volume = []; energy = []
    for i in file[line_from:line_to]:
        scaling_factor.append(float(i.split()[0]))
        volume.append(float(i.split()[1]))
        energy.append(float(i.split()[2]))

    (p, r_squared, volume_fit, energy_fit) = polyfit(volume, energy, 2)
    scaling_factor_fit = np.linspace(scaling_factor[0], scaling_factor[-1], 100)
    plt.plot(volume, energy, '.', volume_fit, energy_fit, '-')
    result_str = "E = %.3f x^2 + (%.3f) x + (%.3f)\nR-squared is %f" % (p[2], p[1], p[0], r_squared)
    plt.text(volume_fit[len(volume_fit)/4], energy_fit[len(energy_fit)/8*1], result_str)
    plt.xlabel(r'Volume ($\AA^{3}$)')
    plt.ylabel('E (eV)')
    print "Fitting result of E-V:", result_str
    print("The scaling factor is %f" % (scaling_factor_fit[energy_fit.argmin()]))
        
elif re.search('.*c[1-9][1-9].*', test_type):     # meaning elastic const.
    delta = []; energy = []
    for i in file[line_from:line_to]:
        delta.append(float(i.split()[0]))
        energy.append(float(i.split()[1]))
        
    (p, r_squared, delta_fit, energy_fit) = polyfit(delta, energy, 2)
    plt.plot(delta, energy, '.', delta_fit, energy_fit, '-')
    result_str = "E = %.3f x^2 + (%.3f) x + (%.3f)\nR-squared is %f" % (p[2], p[1], p[0], r_squared)
    plt.text(delta_fit[len(delta_fit)/4], energy_fit[len(energy_fit)/8*1], result_str)
    plt.xlabel('Delta (ratio)')
    plt.ylabel('E (eV)')
    print "Fitting result:", result_str
    
#    a, c, r_squared, y_fit = centralfit(delta, energy)
#    plt.plot(delta, y_fit, '-', label='E = a*x^2 + c')
#    result_str2 = "E = %.3f x^2 + (%.3f)\nR-squared is %f" % (a, c, r_squared)
#    plt.legend()
#    print "Fitting result:", result_str2

plt.savefig('fit_curve.png')
#plt.show()
