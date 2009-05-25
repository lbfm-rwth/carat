#!/bin/bash
REV_HELP="/home/tilman/c/carat/reverse_help/reverse_help"
REV_HELP2="/home/tilman/c/carat/reverse_help/reverse_help2"

ls | grep lattice > tt

GRPS="`cat tt | sed "s/   / /g" | sed "s/lattice_//g"`"
rm tt
rm pp

# make an intermediate file which does only contain family
# symbols and stuff
for x in $GRPS ; do
   echo x $x 1st

   $REV_HELP $x >> pp

done


for x in $GRPS ; do
   echo x $x
   y=`echo $x | sed "s/_/ /g"`
   echo y $y

   grep "$y" pp > tt
   cat tt | wc -l > reverse_$x

   cat reverse_$x tt | $REV_HELP2 > hh
   mv hh tt

   cat tt | sed "s/inter_//" | sed "s/_/ /g" >> reverse_$x
   rm tt
done

rm pp
