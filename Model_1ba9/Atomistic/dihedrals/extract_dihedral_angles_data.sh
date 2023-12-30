#! /bin/bash
#
for angleFile in *.xvg; do
	gmx analyze -f $angleFile 2> /dev/null \
		| sed -n 7p | awk '{print $2 " " $3}'\
		>> dihedral_mean_stdev.dat
done
