#!/usr/local/bin/bash

ls | grep reverse_ > tt

GR="`cat tt | sed "s/   / /g" | sed "s/reverse_//g"`"

rm $1

for x in $GR ; do
   echo x $x
   y=`echo $x | sed "s/_/ /g"`
   echo y $y

   head -n 1 -l reverse_$x > tt

   sed "s/^1$/xxxxxxx/" tt > pp

   grep xxxx pp > tt

   if [ -s tt ] ; then
      echo $y >> $1
   fi

done

rm -f tt pp
