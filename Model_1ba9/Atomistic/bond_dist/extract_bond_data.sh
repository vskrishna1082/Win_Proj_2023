#! /bin/bash
#
for bondFile in *.xvg; do
	gmx analyze -f $bondFile 2> /dev/null \
		| sed -n 7p | awk '{print $2 " " $3}'\
		>> bond_mean_stdev.dat
done
