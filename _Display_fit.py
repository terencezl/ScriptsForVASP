#!/usr/bin/python
# _Display_fit.py TestType LineFrom LineTo
import sys
import numpy as np
import matplotlib.pyplot as plt

f = open(sys.argv[1]+'_output.txt','rU')
file = f.readlines()
f.close()

line_from = int(sys.argv[2]); line_to = int(sys.argv[3])
x = []; y = []

if sys.argv[1] == 'entest' or sys.argv[1] == 'kptest':
    for i in file[line_from:line_to]:
        x.append(float(i.split()[0]))
        y.append(float(i.split()[3]))
    plt.plot(x, y, '.')
    plt.xlabel('step')
    plt.ylabel('E (eV)')
    
elif sys.argv[1] == 'lctest' or sys.argv[1].isalpha() == False: # meanings elastic const.
    for i in file[line_from:line_to]:
        x.append(float(i.split()[0]))
        y.append(float(i.split()[1]))

    coefficients = np.polyfit(x, y, 2)
    p = np.poly1d(coefficients)
    xp = np.linspace(sorted(x)[0], sorted(x)[-1], 100)
    yp = p(xp)
    plt.plot(x, y, '.', xp, yp, '-')
    
    result_str = "Fitting result: E = %.3f x^2 + (%.3f) x + (%.3f)" % (p[2], p[1], p[0])
    plt.text(xp[len(xp)/4], yp[len(yp)/8*1], result_str)
    print(result_str)
    plt.xlabel('step')
    plt.ylabel('E (eV)')

plt.savefig('fit_curve.png')
#plt.show()
