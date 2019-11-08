#!/bin/sh
set -e

echo "### Test Ex10-1"
../bin/Bravais_catalog << EOF
6-1
y
stdout
s
1
1
n
EOF
echo "### Ex10-1 return code $?"

echo "### Test Ex10-2"
../bin/Elt Ex10_m
echo "### Ex10-2 return code $?"

echo "### Test Ex10-3"
../bin/Orbit -g -S Ex10_m Ex10_B
echo "### Ex10-3 return code $?"

echo "### Test Ex10-4"
../bin/Order Ex10_mB2
echo "### Ex10-4 return code $?"
