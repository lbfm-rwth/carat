#!/usr/local/bin/bash

ls > tt
grep lattice tt > pp
mv pp tt

GROUPS="`cat tt | sed "s/   / /g" | sed "s/lattice_//g"`"

for x in $GROUPS ; do
   echo x $x
   y=`echo $x | sed "s/_/ /g"`
   echo y $y
   rm tt

   # check whether the first group in this file has the same symbol as
   # represented by the filename
   head -l -n3 lattice_$x | grep -l "$y" > tt
   if [ -s tt ] ; then
     echo OK
   else
     echo NOT_OK symbol: $y
   fi

   # check whether the first group appears exactly once in the file
   head -l -n4 lattice_$x | sed "s/^1 $/xxxxxxx/g" | grep -l "xxxxx" > tt
   if [ -s tt ] ; then
     echo OK
   else
     echo NOT_OK symbol: $y
   fi

done

rm -f tt
