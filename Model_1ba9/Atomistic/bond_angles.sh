#! /bin/bash
## Obtai the angles data between set of 3 residues

index=0
firstGroup=
midGroup=
grep -oh "\w*Protein_\w*" index.ndx | while read -r group; do
   if [[ $index == 0 ]]
   then
	   ((index++));
	   firstGroup=$group;
	   continue
   fi
   if [[ $index == 1 ]]
   then
	   ((index++));
	   midGroup=$group;
	   continue
   fi
   gmx gangle -g1 angle -group1 "cog of group $firstGroup plus cog of group $midGroup plus cog of group $group" -n index.ndx -s md_0_1.tpr -f md_0_1_noPBC.xtc -oav "angles/angle$(printf %03d $((index - 1)))";
   firstGroup=$midGroup;
   midGroup=$group;
   ((index++))
done

