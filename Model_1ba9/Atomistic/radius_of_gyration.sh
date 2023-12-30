#! /bin/bash
# Obtains the radius of gyration of each residue

# get all residue groups
grep -oh "\w*Protein_\w*" index.ndx | while read -r group; do
    echo $group | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o "rg/$group.xvg" -n index.ndx
done
