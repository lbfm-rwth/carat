
echo "### Test Ex6-1"
../bin/carat/Bravais_catalog << EOF
1;1;1;1
y
stdout
a
EOF
echo "### Ex6-1 return code $?"

echo "### Test Ex6-2"
../bin/carat/Bravais_catalog << EOF
4-1
y
stdout
a
EOF
echo "### Ex6-2 return code $?"

echo "### Test Ex6-3"
../bin/carat/Bravais_inclusions Ex6_B
echo "### Ex6-3 return code $?"

echo "### Test Ex6-4"
../bin/carat/Tr_bravais Ex6_b3
echo "### Ex6-4 return code $?"

echo "### Test Ex6-5"
../bin/carat/Bravais_type Ex6_b3t
echo "### Ex6-5 return code $?"

echo "### Test Ex6-6"
../bin/carat/Tr_bravais Ex6_b4
echo "### Ex6-6 return code $?"

echo "### Test Ex6-7"
../bin/carat/Bravais_type Ex6_b4t
echo "### Ex6-7 return code $?"
