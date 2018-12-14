
../bin/carat/Symbol Ex2_G

../bin/carat/Bravais_catalog << EOF
3;2-2;1
n
EOF

../bin/carat/Order -o Ex2_G

../bin/carat/QtoZ -D Ex2_Go

for f in Ex2_Go.*; do ../bin/carat/Bravais_type $f; done
rm -f Ex2_Go.*
