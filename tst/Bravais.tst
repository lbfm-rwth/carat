# This is a CARAT test file. It must have extension `.tst`. 
# Comments and blank lines are ignored. For a test file named 
# `testfile.tst`, the correct output of `bash testfile.tst`
# must be placed in the file `testfile.out`.

# Testing Form_space
echo "### Test Bravais-1"
../bin/carat/Form_space grp1.dat
echo "### Bravais-1 return code $?"

# Testing Bravais_grp
echo "### Test Bravais-2"
../bin/carat/Bravais_grp grp1.dat
echo "### Bravais-2 return code $?"

