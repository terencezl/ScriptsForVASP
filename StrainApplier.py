#!/usr/bin/env python
# Used by SequenceTest.sh for altering the real space basis_vectors in POSCAR

import sys
import numpy as np
import argparse


def main(arguments='-h'):
    """
    The main function that does the strain application.

    :param arguments: string
        Command line arguments and options. Typically ' '.join(sys.argv[1:]).

    Usages
    ------
    main(), main('-h') : as if ``StrainApplier.py -h`` is executed from command line.
    main(string) : as if ``StrainApplier.py content_of_string`` is executed from command line.
    """
    # parse arguments
    arguments = arguments.split()
    parser = argparse.ArgumentParser(description="""Apply strain to the cell basis vectors in POSCAR.""")
    parser.add_argument('test_type', metavar='TEST_TYPE', help="the elastic constant combination to select, like c11+2c12")
    parser.add_argument('cryst_sys', metavar='CRYST_SYS', help="the crystallographic system to select, like cubic")
    parser.add_argument('delta', type=float, help="the delta value for the strain")
    parser.add_argument('-i', '--input', metavar='POSCAR', default='POSCAR', help="the input POSCAR file name")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o', '--output', metavar='POSCAR-strained', default='default', help="the output file name")
    group.add_argument('-p', '--inplace', action='store_true',
                                    help="directly change POSCAR in place without creating a new file")
    args = parser.parse_args(arguments)
    delta = args.delta
    if args.output == 'default':
        args.output = args.input + '-strained'
    if args.inplace:
        args.output = args.input

    # open the POSCAR file
    with open(args.input, 'r') as f:
        POSCAR = f.readlines()
    # obtain the basis vector 2D array
    basis_vectors = np.zeros((3, 3))
    for i, line in enumerate(POSCAR[2:5]):
        basis_vectors[i] = line.split()

    if args.cryst_sys == 'cubic':
        if args.test_type == "c11+2c12":
            transformation_matrix = np.array([[1 + delta, 0, 0],
                                              [0, 1 + delta, 0],
                                              [0, 0, 1 + delta]])
        elif args.test_type == "c11-c12":
            transformation_matrix = np.array([[1 + delta, 0, 0],
                                              [0, 1 - delta, 0],
                                              [0, 0, 1 + delta ** 2 / (1 - delta ** 2)]])

        elif args.test_type == "c44":
            transformation_matrix = np.array([[1, delta / 2., 0],
                                              [delta / 2., 1, 0],
                                              [0, 0, 1 + delta ** 2 / (4 - delta ** 2)]])

    elif args.cryst_sys == 'tetragonal':
        if args.test_type == "c11":
            transformation_matrix = np.array([[1 + delta, 0, 0],
                                              [0, 1, 0],
                                              [0, 0, 1]])
        elif args.test_type == "c33":
            transformation_matrix = np.array([[1, 0, 0],
                                              [0, 1, 0],
                                              [0, 0, 1 + delta]])
        elif args.test_type == "c44":
            transformation_matrix = np.array([[1, 0, 0],
                                              [0, 1, delta],
                                              [0, delta, 1]])
        elif args.test_type == "5c11-4c12-2c13+c33":
            transformation_matrix = np.array([[1 + 2 * delta, 0, 0],
                                              [0, 1 - delta, 0],
                                              [0, 0, 1 - delta]])
        elif args.test_type == "c11+c12-4c13+2c33":
            transformation_matrix = np.array([[1 - delta, 0, 0],
                                              [0, 1 - delta, 0],
                                              [0, 0, 1 + 2 * delta]])
        elif args.test_type == "c11+c12-4c13+2c33+2c66":
            transformation_matrix = np.array([[1 + delta, delta, 0],
                                              [delta, 1 + delta, 0],
                                              [0, 0, 1 - 2 * delta]])

    elif args.cryst_sys == 'orthorhombic':
        if args.test_type == "c11":
            transformation_matrix = np.array([[1 + delta, 0, 0],
                                              [0, 1, 0],
                                              [0, 0, 1]])
        elif args.test_type == "c22":
            transformation_matrix = np.array([[1, 0, 0],
                                              [0, 1 + delta, 0],
                                              [0, 0, 1]])
        elif args.test_type == "c33":
            transformation_matrix = np.array([[1, 0, 0],
                                              [0, 1, 0],
                                              [0, 0, 1 + delta]])
        elif args.test_type == "c44":
            transformation_matrix = np.array([[1, 0, 0],
                                              [0, 1, delta / 2],
                                              [0, delta / 2, 1]])
        elif args.test_type == "c55":
            transformation_matrix = np.array([[1, 0, delta / 2],
                                              [0, 1, 0],
                                              [delta / 2, 0, 1]])
        elif args.test_type == "c66":
            transformation_matrix = np.array([[1, delta / 2, 0],
                                              [delta / 2, 1, 0],
                                              [0, 0, 1]])
        elif args.test_type == "4c11-4c12-4c13+c22+2c23+c33":
            transformation_matrix = np.array([[1 + 2 * delta, 0, 0],
                                              [0, 1 - delta, 0],
                                              [0, 0, 1 - delta]])
        elif args.test_type == "c11-4c12+2c13+4c22-4c23+c33":
            transformation_matrix = np.array([[1 - delta, 0, 0],
                                              [0, 1 + 2 * delta, 0],
                                              [0, 0, 1 - delta]])
        elif args.test_type == "c11+2c12-4c13+c22-4c23+4c33":
            transformation_matrix = np.array([[1 - delta, delta, 0],
                                              [delta, 1 - delta, 0],
                                              [0, 0, 1 + 2 * delta]])

    # apply the specified strain
    basis_vectors_new = np.dot(basis_vectors, transformation_matrix)
    # create the string list form, ready to write to file
    basis_vectors_new_str = [''] * len(basis_vectors_new)
    for i in range(3):
        basis_vectors_new_str[i] = '{0[0]:11.8f}  {0[1]:11.8f}  {0[2]:11.8f}\n'.format(basis_vectors_new[i])
    # replace the POSCAR variable content
    for i, line_num in enumerate(range(2, 5)):
        POSCAR[line_num] = basis_vectors_new_str[i]
    # write to it
    with open(args.output, 'w') as f:
        f.writelines(POSCAR)


if __name__ == '__main__':
    main(' '.join(i.replace(' ', '') for i in sys.argv[1:]))