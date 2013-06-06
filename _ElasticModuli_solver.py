#!/usr/local/python/2.7.1/bin/python
# Used by ElasticSolve.sh to solve linear equations to get elastic consts
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# _ElasticModuli_solver.py cryst_sys volumn_of_primitive_cell input_data

import sys
import numpy as np

cryst_sys = sys.argv[1]
Vpcell = float(sys.argv[2])
N_atoms = np.array(eval(sys.argv[4]))
econst = np.array(eval(sys.argv[4]))

total_atoms = np.sum(N_atoms)
N_A = 6.0221413e+23
k_B = 1.3806488e-23
h = 6.62606957e-34

if cryst_sys == 'cubic':
    [c11, c12, c44] = econst
    B = (c11 + 2*c12)/3.
    G_V = ((c11 - c12) + 3*c44)/5.
    G_R = 5 * (c11 - c12) * c44 / (4*c44 + 3*c11 - 3*c12)
    G = (G_V + G_R)/2.
    E = 9*B*G/(3*B + G)
    nu = (3*B - 2*G)/(2 * (3*B + G))
    k = G/B
    cauchy_pressure = c12 - c44
    H_V = (1 - 2*nu) * E/(6 * (1 + nu))
    row = (N_atoms[0] * 44.96 + N_atoms[1] * 14.01) * 10**27 / Vpcell / N_A
    v_l = ((B + 4/3. * G)/row)**0.5
    v_t = (G/row)**0.5
    v_m = (1/3. * (2/v_t**3 + 1/v_l**3))**(-1/3.)
#    row = (N_atoms[0] * atomic_mass_M + N_atoms[1] * atomic_mass_N) / Vpcell / N_A
    theta_D = h/k_B * (3*total_atoms/4/np.pi/Vpcell)**(1/3.) * v_m
    print("B = %f\nG_V = %f\nG_R = %f\nG = %f\nE = %f\nnu = %f\nk = %f\ncauchy_pressure = %f\nH_V = %f\ntheta_D = %f" % (B, G_V, G_R, G, E, nu, k, cauchy_pressure, H_V, theta_D))

if cryst_sys == 'tetragonal':
    pass

if cryst_sys == 'orthorhombic':
    pass
    