
echo "### Test Ex12-1"
../bin/carat/Q_catalog << EOF
s abb min.108
w min.108
q
EOF
echo "### Ex12-1 return code $?"
cat min.108

echo "### Test Ex12-2"
../bin/carat/QtoZ -D Ex12_min.108
echo "### Ex12-2 return code $?"

for x in Ex12_min.108.?.? ; do
  echo "### Test Ex12-3-$x"
  ../bin/carat/Extensions $x > ex.$x
  echo "### Ex12-3-$x return code $?"
done
cat ex.*

echo "### Test Ex12-4"
../bin/carat/Presentation Ex12_min.108 > pres
echo "### Ex12-4 return code $?"
cat pres
for x in Ex12_min.108.?.? ; do
  echo "### Test Ex12-5-$x"
  ../bin/carat/Extensions pres $x > ex2.$x
  echo "### Ex12-5-$x return code $?"
done
cat ex2.*

for x in Ex12_min.108.?.?.? ; do
  echo "### Test Ex12-6-$x"
  ../bin/carat/Torsionfree $x
  echo "### Ex12-6-$x return code $?"
  echo
done

rm -f min.108 ex.* ex2.* pres Ex12_min.108.?.?
