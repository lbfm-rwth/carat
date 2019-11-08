#!/bin/sh
set -e

echo "### Test Ex7-1"
../bin/Normalizer Ex7_g
echo "### Ex7-1 return code $?"

echo "### Test Ex7-2"
../bin/Is_finite Ex7_nh
echo "### Ex7-2 return code $?"

echo "### Test Ex7-3"
../bin/Symbol Ex7_g
echo "### Ex7-3 return code $?"
