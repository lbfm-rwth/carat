
echo "### Test Ex14-1"
../bin/carat/KSubgroups -a -t -n Ex14_R 2
echo "### Ex14-1 return code $?"

echo "### Test Ex14-2"
../bin/carat/KSupergroups -a -t -n Ex14_R 2
echo "### Ex14-2 return code $?"

echo "### Test Ex14-3"
../bin/carat/TSubgroups -a -t Ex14_R
echo "### Ex14-3 return code $?"

# this one segfaults:
# echo "### Test Ex14-4"
# ../bin/carat/TSupergroups -t Ex14_R
# echo "### Ex14-4 return code $?"

# who creates this?
rm -f ZZ.tmp
