
../bin/carat/Bravais_catalog << EOF
1;1;1;1
y
stdout
a
EOF

../bin/carat/Bravais_catalog << EOF
4-1
y
stdout
a
EOF

../bin/carat/Bravais_inclusions Ex6_B
grep 1\;1\;1\;1 Ex6_Bin

../bin/carat/Tr_bravais Ex6_b3
../bin/carat/Bravais_type Ex6_b3t
../bin/carat/Tr_bravais Ex6_b4
../bin/carat/Bravais_type Ex6_b4t
