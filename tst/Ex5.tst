#!/bin/sh
set -e

echo "### Test Ex5-1"
../bin/Extract Ex5_g
echo "### Ex5-1 return code $?"

echo "### Test Ex5-2"
../bin/Is_finite Ex5_gp
echo "### Ex5-2 return code $?"

echo "### Test Ex5-3"
../bin/Presentation Ex5_gp
echo "### Ex5-3 return code $?"

echo "### Test Ex5-4"
../bin/Standard_affine_form Ex5_g Ex5_gpp
echo "### Ex5-4 return code $?"

echo "### Test Ex5-5"
../bin/Extract Ex5_C
echo "### Ex5-5 return code $?"

echo "### Test Ex5-6"
../bin/Order Ex5_P
echo "### Ex5-6 return code $?"

echo "### Test Ex5-7"
../bin/Presentation Ex5_P
echo "### Ex5-7 return code $?"

echo "### Test Ex5-8"
../bin/Standard_affine_form Ex5_C Ex5_Ppres
echo "### Ex5-8 return code $?"

echo "### Test Ex5-9"
../bin/Name -o Ex5_S
echo "### Ex5-9 return code $?"
