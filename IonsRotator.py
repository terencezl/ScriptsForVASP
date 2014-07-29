#!/usr/bin/env python
# IonsRotator.py -l [23,24,25,26,27,28] -c [0,0.5,0.75] -u [np.sqrt(2),0,1] -a 60 -i POSCAR -o POSCAR-rotated

import sys
import numpy as np
import argparse


def main(arguments='-h'):
    """
    The main function that does the ion rotation.

    :param arguments: string
        Command line arguments and options. Typically ' '.join(sys.argv[1:]).

    Usages
    ------
    main(), main('-h') : as if ``IonsRotator.py -h`` is executed from command line.
    main(string) : as if ``IonsRotator.py content_of_string`` is executed from command line.
    """
    arguments = arguments.split()
    parser = argparse.ArgumentParser(description="""Rotate the ions in certain lines of POSCAR
                            around a certain point, along a certain direction, to a certain angle.""")
    parser.add_argument('-l', '--lines', type=eval, required=True, help="""the line numbers in
                            POSCAR to operate on in the form of '[11,12,13,14,...]'. Count from 1.""")
    parser.add_argument('-c', '--center', type=eval, required=True, help="""the rotation center's
                            direct coordinate in the form of '(a_1,a_2,a_3)]'""")
    parser.add_argument('-u', '--direction', type=eval, required=True, help="""the rotation direction's
                            cartesian coordinate in the form of '(u_x,u_y,u_z)'. No need to be normalized.
                            Values can be numpy functions like np.sqrt(2).""")
    parser.add_argument('-a', '--angle', type=float, required=True, help="""the rotation angle in degree,
                            counter-clockwise""")
    parser.add_argument('-i', '--input', metavar='POSCAR', default='POSCAR', help="the input POSCAR file name")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o', '--output', metavar='POSCAR-rotated', default='default', help="the output file name")
    group.add_argument('-p', '--inplace', action='store_true',
                                    help="directly change POSCAR in place without creating a new file")
    args = parser.parse_args(arguments)
    if args.output == 'default':
        args.output = args.input + '-rotated'
    if args.inplace:
        args.output = args.input

    # np.set_printoptions(suppress=True)
    ions_line_number = np.array(args.lines) - 1  # machine counts from 0

    # define the rotation matrix A. Here we use the active rotation scheme, i.e. rotating the object
    [u_x, u_y, u_z] = args.direction / np.linalg.norm(args.direction)
    angle = args.angle / 180. * np.pi
    A = np.eye(3) * np.cos(angle) + np.array([[0, -u_z, u_y], [u_z, 0, -u_x], [-u_y, u_x, 0]]) \
                                    * np.sin(angle) + np.array([[u_x ** 2, u_x * u_y, u_x * u_z],
                                                                [u_x * u_y, u_y ** 2, u_y * u_z],
                                                                [u_x * u_z, u_y * u_z, u_z ** 2]]) * (1 - np.cos(angle))

    # open the POSCAR that has ions' position part
    with open(args.input, 'r') as f:
        POSCAR = f.readlines()

    # obtain the basis vectors from POSCAR
    basis_vectors = np.zeros((3, 3))
    for i, line in enumerate(POSCAR[2:5]):
        basis_vectors[i] = line.split()

    # convert the string form of ions' position into a numpy numeral array and subtract the translational vector
    ions_position = np.zeros((len(ions_line_number), 3))
    for i, line in enumerate(ions_line_number):
        ions_position[i] = POSCAR[line].split()[:3]
    ions_position -= args.center
    # transform the coordinates from direct to cartesian
    ions_position = np.dot(ions_position, basis_vectors)
    # the matrix multiplication
    ions_position_new = np.dot(ions_position, A.transpose())
    # transform the coordinates from cartesian back to direct
    ions_position_new = np.dot(ions_position_new, np.linalg.inv(basis_vectors))

    # ions_position_new = ions_position_new / np.cos(angle)

    # add the translational vector
    ions_position_new = ions_position_new + args.center

    ions_position_new_str = [''] * len(ions_position_new)
    for i in range(len(ions_position_new)):
        ions_position_new_str[i] = '{0[0]:11.8f}  {0[1]:11.8f}  {0[2]:11.8f}\n'.format(ions_position_new[i])

    for i, line in enumerate(ions_line_number):
        POSCAR[line] = ions_position_new_str[i]

    with open(args.output, 'w') as f:
        f.writelines(POSCAR)


if __name__ == '__main__':
    main(' '.join(i.replace(' ', '') for i in sys.argv[1:]))