#!/bin/sh
set -e

echo "### Test Ex14-1"
../bin/KSubgroups -a -t -n Ex14_R 2
echo "### Ex14-1 return code $?"

echo "### Test Ex14-2"
../bin/KSupergroups -a -t -n Ex14_R 2
echo "### Ex14-2 return code $?"

echo "### Test Ex14-3"
../bin/TSubgroups -a -t Ex14_R
echo "### Ex14-3 return code $?"

echo "### Test Ex14-4"
../bin/TSupergroups -t Ex14_R
echo "### Ex14-4 return code $?"

# who creates this?
rm -f ZZ.tmp
