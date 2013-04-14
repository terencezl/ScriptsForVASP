mkdir $1 2> /dev/null
cd $1
cp ../INPUT/INCAR .
cp ../INPUT/POSCAR .
cp ../INPUT/POTCAR .
cp ../INPUT/KPOINTS .
cp ../INPUT/qsub.parallel .
qname=${PWD//\//_}                                                  # edit the name of the command so that it looks like terencelz_GaN_MN_lctest_140
qname=${qname##*utl0268_}
sed -i s/@N@/$qname/g qsub.parallel                                 # replace the name in the command file. Arg reserved
sed -i s%@R@%$PWD%g qsub.parallel                                   # replace the trial subfolder in the command file. Arg reserved