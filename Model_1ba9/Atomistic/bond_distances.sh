#! /bin/bash
# Obtain a list of bonds between residues

index=0
oldGroup=
grep -oh "\w*Protein_\w*" index.ndx | while read -r group; do
    if [[ $index == 0 ]] # bond starts with 1-2
    then
	    ((index++));
	    oldGroup=$group;
	    continue
    fi
    gmx distance -s md_0_1.tpr -f md_0_1_noPBC.xtc -n index.ndx -oav "bond_dist/bond$(printf %03d $index)" -select "cog of group $group plus cog of group $oldGroup";
    ((index++));
    oldGroup=$group;
done
