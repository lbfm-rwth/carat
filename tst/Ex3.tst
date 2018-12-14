
../bin/carat/Bravais_catalog << EOF
1,1,1,1,1
y
stdout
a
EOF

../bin/carat/Bravais_inclusions Ex3_11111 -S > all

../bin/carat/Bravais_catalog << EOF
5-1
y
stdout
a
EOF

../bin/carat/Bravais_catalog << EOF
5-2
y
stdout
a
EOF

../bin/carat/Bravais_inclusions Ex3_51a  > notmax
../bin/carat/Bravais_inclusions Ex3_51b >> notmax
../bin/carat/Bravais_inclusions Ex3_51c >> notmax
../bin/carat/Bravais_inclusions Ex3_52a >> notmax
../bin/carat/Bravais_inclusions Ex3_52b >> notmax
../bin/carat/Bravais_inclusions Ex3_52c >> notmax
../bin/carat/Bravais_inclusions Ex3_52d >> notmax
grep Symbol notmax | sort -u

sort all > allsort
diff allsort Ex3_MAX | grep Symbol

rm -f all allsort notmax
