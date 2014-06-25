#!/usr/bin/env python
# _Display_fit.py test_type line_begin line_to
import sys
import numpy as np
import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpltools.style
mpltools.style.use('ggplot')
from scipy.optimize import curve_fit
import re


def polyfit(x, y, degree):
    popts = np.polyfit(x, y, degree, full=True)
    p = np.poly1d(popts[0])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 1000)
    y_fit = p(x_fit)
    y_avg = np.sum(y) / len(y)
    ss_total = np.sum((y - y_avg) ** 2)
    r_squared = 1 - popts[1] / ss_total
    return p, r_squared[0], x_fit, y_fit


#def murnaghan_eqn(V, V0, B0, B0_prime, E0):
#    return (E0 + B0/B0_prime * V * (1 + (V0/V)**B0_prime/(B0_prime - 1)) - V0 * B0/(B0_prime -1))

def murnaghan_eqn(V, V0, B0, B0_prime, E0):
    return (E0 + 9 * V0 * B0 / 16. * (
    ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime + ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.))))


def murnaghan_fit(x, y):
    coeffs, pcov = curve_fit(murnaghan_eqn, x, y, [np.average(x), 2.5, 4, np.average(y)])
    x_fit = np.linspace(sorted(x)[0], sorted(x)[-1], 1000)
    y_fit = murnaghan_eqn(x_fit, *coeffs)
    y_fit_eqlen = murnaghan_eqn(x, *coeffs)
    ss_resid = np.sum((y_fit_eqlen - y) ** 2)
    y_avg = np.sum(y) / len(y)
    ss_total = np.sum((y - y_avg) ** 2)
    r_squared = 1 - ss_resid / ss_total
    return coeffs, r_squared, x_fit, y_fit


test_type = sys.argv[1]
line_begin = int(sys.argv[2])
data_line_count = int(sys.argv[3])
with open(test_type + '_output.txt', 'rU') as f:
    test_output = f.readlines()
for i, line in enumerate(test_output):
    test_output[i] = line.split()

if test_type == 'entest':
    col_names = ['ENCUT', 'E_LC1(eV)', 'dE_LC1 (eV)', 'E_LC2(eV)', 'dE_LC2 (eV)', 'DE=E_LC2-E_LC1(eV)', 'dDE(eV)']
    data = np.zeros((data_line_count, 7))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    plt.plot(data[:, 0], data[:, 2], 'o', label=col_names[2])
    plt.plot(data[:, 0], data[:, 4], 's', label=col_names[4])
    plt.plot(data[:, 0], data[:, 6], '^', label=col_names[6])
    plt.axhline(0.001, ls=':', c='k')
    plt.grid(True)
    plt.xlabel('ENCUT (eV)')
    plt.ylabel('Energy diff (eV)')

elif test_type == 'kptest':
    col_names = ['nKP', 'E_LC1 (eV)', 'dE_LC1 (eV)', 'E_LC2 (eV)', 'dE_LC2 (eV)', 'DE=E_LC2-E_LC1 (eV)', 'dDE (eV)']
    data = np.zeros((data_line_count, 7))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    plt.plot(data[:, 0], data[:, 2], 'o', label=col_names[2])
    plt.plot(data[:, 0], data[:, 4], 's', label=col_names[4])
    plt.plot(data[:, 0], data[:, 6], '^', label=col_names[6])
    plt.axhline(0.001, ls=':', c='k')
    plt.grid(True)
    plt.xlabel('nKP')
    plt.ylabel('Energy diff (eV)')

elif test_type == 'lctest':
    col_names = ['Scaling factor (Ang)', 'Volume (Ang^3)', 'E (eV)']
    data = np.zeros((data_line_count, 3))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    V_a_conversion_multiplier = data[0, 0] ** 3 / data[0, 1]
    plt.plot(data[:, 1], data[:, 2], 'o')

    # fitting the Birch-Murnaghan equation of state
    (coeffs_vol_M, r_squared_M, volume_fit_M, energy_fit_M) = murnaghan_fit(data[:, 1], data[:, 2])
    scaling_factor_eqlbrm = (coeffs_vol_M[0] * V_a_conversion_multiplier) ** (1 / 3.)

    # plotting the Birch-Murnaghan equation of state
    plt.plot(volume_fit_M, energy_fit_M, '-', label="B-M eqn of state")
    result_str = "R-squared is %f" % r_squared_M
    plt.text(volume_fit_M[len(volume_fit_M) / 4], energy_fit_M[60], result_str)

    # standrad output, directed to files by the bash script calling this python script
    print "%s" % result_str
    print("Equilibrium scaling factor is {0}".format(scaling_factor_eqlbrm))
    if scaling_factor_eqlbrm <= data[0, 0] or scaling_factor_eqlbrm >= data[-1, 0]:
        print("!Equilibrium point is out of the considered range!")
    else:
        print("V0 = %f\nB0 = %f\nB0' = %f" % (coeffs_vol_M[0], coeffs_vol_M[1] * 160.2, coeffs_vol_M[2]))
        print("Total energy is {0}".format(coeffs_vol_M[3]))
        np.savetxt('eosfit_data.dat', np.column_stack((volume_fit_M, energy_fit_M)),
                   '%15.6E', header='Volume (Ang^3) E (eV) ')

    plt.xlabel(r'Volume ($\AA^{3}$)')
    plt.ylabel('E (eV)')
    plt.legend(loc=0)
    np.savetxt('orig_data.dat', data, '%15.6E', header=' '.join(col_names))

elif test_type == 'rttest':
    col_names = ['Ratio', 'Volume (Ang^3)', 'E (eV)']
    data = np.zeros((data_line_count, 3))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    plt.plot(data[:, 0], data[:, 2], 'o')
    (p_ratio, r_squared, ratio_fit, energy_fit) = polyfit(data[:, 0], data[:, 2], 3)
    ratio_eqlbrm = ratio_fit[energy_fit.argmin()]
    plt.plot(ratio_fit, energy_fit, '-')

    print("R-squared is {0}\nEquilibrium ratio is {1}\nThe polynomial is\n{2}".format(r_squared, ratio_eqlbrm, p_ratio))
    print("Minimal total energy is %f" % energy_fit.min())
    np.savetxt('polyfit_data.dat', np.column_stack((ratio_fit, energy_fit)),
               '%15.6E', header='Ratio E (eV) ')
    plt.xlabel('Raito')
    plt.ylabel('E (eV)')
    np.savetxt('orig_data.dat', data, '%15.6E', header=' '.join(col_names))

elif test_type == 'agltest':
    col_names = ['Angle (degree)', 'E (eV)']
    data = np.zeros((data_line_count, 2))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    plt.plot(data[:, 0], data[:, 1], 'o')
    plt.xlabel('Angle (degree)')
    plt.ylabel('E (eV)')
    np.savetxt('orig_data.dat', data, '%15.6E', header=' '.join(col_names))

elif re.search('.*c[1-9][1-9].*', test_type) or re.search('A.*', test_type):
    col_names = ['Delta', 'E (eV)']
    data = np.zeros((data_line_count, 2))
    for i, row in test_output[line_begin:line_begin + data_line_count]:
        data[i] = row

    if re.search('.*c[1-9][1-9].*', test_type):
        (p, r_squared, delta_fit, energy_fit) = polyfit(data[:, 0], data[:, 1], 2)
        plt.plot(data[:, 0], data[:, 1], 'o', delta_fit, energy_fit, '-')
        result_str = "E = %f x^2 + (%f) x + (%f)\nR-squared is %f" % (p[2], p[1], p[0], r_squared)
    elif re.search('A.*', test_type):
        (p, r_squared, delta_fit, energy_fit) = polyfit(data[:, 0], data[:, 1], 3)
        plt.plot(data[:, 0], data[:, 1], 'o', delta_fit, energy_fit, '-')
        result_str = "E = %f x^3 + (%f) x^2 + (%f) x + (%f)\nR-squared is %f" % (p[3], p[2], p[1], p[0], r_squared)

    plt.text(delta_fit[len(delta_fit) / 4], energy_fit[6], result_str)
    plt.xlabel('Delta (ratio)')
    plt.ylabel('E (eV)')
    print "Fitting result:", result_str

plt.savefig('fit_curve.png')