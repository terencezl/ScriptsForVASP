#!/usr/bin/python
import sys
#import numpy as np
import matplotlib.pyplot as plt

N_step = 301
if len(sys.argv) != 2:
    N_step = int(sys.argv[-1])

f = open('DOSCAR','rU')
list = [];dos = [];dosx = [];dosy_tot = []; dosy_int = []; dosy_s = []; dosy_p =[]; dosy_d = []
for line in f:
    list.append(line[0:-1])
f.close()

# the total DOS and the integrated DOS
fenergy = float(list[5].split()[3])
for i in range( 6, 6 + N_step ):
    dos.append(list[i].split())
for i in dos:
    dosx.append(float(i[0]) - fenergy)
    dosy_tot.append(float(i[1]))
    dosy_int.append(float(i[2]))
#plt.plot(dosx,dosy_tot, label="total")
plt.plot(dosx,dosy_int, label="integral")

# the projected DOS
for n in range(1, int(sys.argv[1])+1):
    dos = []; dosx = []; dosy_s = [];dosy_p = []; dosy_d = []
    for i in range ( 6 + (N_step+1)*n , 6 + ((N_step+1)*(n+1)-1) ):
        dos.append(list[i].split())
    for i in dos:
        dosx.append(float(i[0]) - fenergy)
        dosy_s.append(float(i[1]))
        dosy_p.append(float(i[2]))
        dosy_d.append(float(i[3]))
    plt.plot(dosx,dosy_s, label=str(n)+"_s")
    plt.plot(dosx,dosy_p, label=str(n)+"_p")
    plt.plot(dosx,dosy_d, label=str(n)+"_d")
    plt.legend()

#plt.savefig(sys.argv[1])
plt.savefig('dos.png')
plt.show()
