echo "### Test Ex10-1"
../bin/carat/Bravais_catalog << EOF
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
../bin/carat/Elt Ex10_m
echo "### Ex10-2 return code $?"

echo "### Test Ex10-3"
../bin/carat/Orbit -g -S Ex10_m Ex10_B
echo "### Ex10-3 return code $?"

echo "### Test Ex10-4"
../bin/carat/Order Ex10_mB2
echo "### Ex10-4 return code $?"
