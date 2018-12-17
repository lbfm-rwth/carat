
echo "### Test Ex2-1"
../bin/carat/Symbol Ex2_G
echo "### Ex2-1 return code $?"

echo "### Test Ex2-2"
../bin/carat/Bravais_catalog << EOF
3;2-2;1
n
EOF
echo "### Ex2-2 return code $?"

echo "### Test Ex2-3"
../bin/carat/Order -o Ex2_G
echo "### Ex2-3 return code $?"

echo "### Test Ex2-4"
../bin/carat/QtoZ -D Ex2_Go
echo "### Ex2-4 return code $?"

for f in Ex2_Go.*; do
  echo "### Test Ex2-5-$f"
  ../bin/carat/Bravais_type $f;
  echo "### Ex2-5-$f return code $?"
done

rm -f Ex2_Go.*
