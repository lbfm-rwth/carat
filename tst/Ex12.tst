
../bin/carat/Q_catalog << EOF
s abb min.108
w min.108
q
EOF
cat min.108

../bin/carat/QtoZ -D Ex12_min.108

for x in Ex12_min.108.?.? ; do
  ../bin/carat/Extensions $x > ex.$x
done
cat ex.*

../bin/carat/Presentation Ex12_min.108 > pres
cat pres
for x in Ex12_min.108.?.? ; do
  ../bin/carat/Extensions pres $x > ex2.$x
done
cat ex2.*

for x in Ex12_min.108.?.?.? ; do
  echo $x
  ../bin/carat/Torsionfree $x
  echo
done

rm -f min.108 ex.* ex2.* pres Ex12_min.108.?.?
