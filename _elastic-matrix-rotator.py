#!/usr/local/bin/python
import numpy as np
np.set_printoptions(suppress=True)

# suffices transformation
def C_to_c(i):
  if i == 0:
    m = (0,0)
  elif i == 1:
    m = (1,1)
  elif i == 2:
    m = (2,2)
  elif i == 3:
    m = (1,2)
  elif i == 4:
    m = (2,0)
  elif i == 5:
    m = (0,1)
  return m

# rotation matrix a[i,j]
a = np.array([[1/np.sqrt(2), 1/np.sqrt(2), 0],
              [-1/np.sqrt(2), 1/np.sqrt(2), 0],
              [0, 0, 1]])

# assigning values to C[i,j]
C = np.zeros((6,6))
C_prime = np.zeros((6,6))
C[0,0] = 384.187668
C[1,1] = C[0,0]
C[2,2] = C[0,0]
C[0,1] = 120.960484
C[1,0] = C[0,1]
C[0,2] = C[0,1]
C[2,0] = C[0,2]
C[1,2] = C[0,1]
C[2,1] = C[0,1]
C[3,3] = 243.561463
C[4,4] = C[3,3]
C[5,5] = C[3,3]

# assigning values automatically to the lower left half
for i in np.arange(6):
  for j in np.arange(i):
    C[i,j] = C[j,i]

# transformation from C[i,j] to c[i,j,k,l]
c = np.zeros((3,3,3,3))
c_prime = np.zeros((3,3,3,3))
#c[0,0,0,0] = C[0,0]
#c[1,1,1,1] = C[1,1]
#c[2,2,2,2] = C[2,2]
#c[0,0,1,1] = C[0,1]
#c[0,0,2,2] = C[0,2]
#c[1,1,2,2] = C[1,2]
#c[1,2,1,2] = c[1,2,2,1] = c[2,1,1,2] = c[2,1,2,1] = C[3,3]
#c[2,0,2,0] = c[2,0,0,2] = c[0,2,2,0] = c[0,2,0,2] = C[4,4]
#c[0,1,0,1] = c[0,1,1,0] = c[1,0,0,1] = c[1,0,1,0] = C[5,5]

for i in np.arange(6):
  for j in np.arange(6):
    m = C_to_c(i)
    n = C_to_c(j)
    c[m + n] = c[m[::-1] + n] = c[m + n[::-1]] = c[m[::-1] + n[::-1]] = C[i,j]

# rotation ransformation from c[i,j,k,l] to c_prime[i,j,k,l]
for i in np.arange(3):
  for j in np.arange(3):
    for k in np.arange(3):
      for l in np.arange(3):
        for m in np.arange(3):
          for n in np.arange(3):
            for o in np.arange(3):
              for p in np.arange(3):
                 c_prime[i,j,k,l] += a[i,m] * a[j,n] * a[k,o] * a[l, p] * c[m,n,o,p]

# transformation from c_prime[i,j,k,l] to C_prime[i,j]
for i in np.arange(6):
  for j in np.arange(6):
    m = C_to_c(i)
    n = C_to_c(j)
    C_prime[i,j] = c_prime[m + n]

# print it out
print('Original elastic matrix:\n{0}\n\nRotated elastic matrix:\n{1}'.format(C, C_prime))
