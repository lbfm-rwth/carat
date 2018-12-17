
echo "### Test Ex1-1"
../bin/carat/Bravais_catalog << EOF
1,1,1,1
y
stdout
a
EOF
echo "### Ex1-1 return code $?"

echo "### Test Ex1-2"
../bin/carat/Bravais_inclusions -S Ex1_1111
echo "### Ex1-2 return code $?"
