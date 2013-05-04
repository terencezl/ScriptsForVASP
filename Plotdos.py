#!/usr/bin/python
import sys
#import numpy as np
import matplotlib.pyplot as plt

f = open('DOSCAR','rU')
list = [];dos = [];dosx = [];dosy = []
for line in f:
    list.append(line[0:-1])
f.close()
fenergy = float(list[5].split()[3])
for i in range(6,507):
    dos.append(list[i].split())
for i in dos:
    dosx.append(float(i[0]) - fenergy)
for i in dos:
    dosy.append(float(i[1]))
plt.plot(dosx,dosy)
#plt.savefig(sys.argv[1])
plt.savefig('dos.png')
plt.show()
