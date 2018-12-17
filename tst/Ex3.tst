
echo "### Test Ex3-1"
../bin/carat/Bravais_catalog << EOF
1,1,1,1,1
y
stdout
a
EOF
echo "### Ex3-1 return code $?"

echo "### Test Ex3-2"
../bin/carat/Bravais_inclusions Ex3_11111 -S
echo "### Ex3-2 return code $?"

echo "### Test Ex3-3"
../bin/carat/Bravais_catalog << EOF
5-1
y
stdout
a
EOF
echo "### Ex3-3 return code $?"

echo "### Test Ex3-4"
../bin/carat/Bravais_catalog << EOF
5-2
y
stdout
a
EOF
echo "### Ex3-4 return code $?"

for f in Ex3_51? ; do
  echo "### Test Ex3-4-$f"
  ../bin/carat/Bravais_inclusions $f
  echo "### Ex3-4-$f return code $?"
done

