
# compared to Ex4_gn, the result has a different set of normalizer generators,
# but the normalizer is the same
echo "### Test Ex4-1"
../bin/carat/Normalizer Ex4_g
echo "### Ex4-1 return code $?"

echo "### Test Ex4-2"
../bin/carat/Vector_systems Ex4_gn
echo "### Ex4-2 return code $?"

echo "### Test Ex4-3"
../bin/carat/Extract -r Ex4_gn Ex4_out
echo "### Ex4-3 return code $?"

echo "### Test Ex4-4"
../bin/carat/Tr_bravais Ex4_gn
echo "### Ex4-4 return code $?"

echo "### Test Ex4-5"
../bin/carat/Sublattices -b Ex4_gn_tr
echo "### Ex4-5 return code $?"

for i in 2 3 4 5 6 ; do
  for j in 2 3 4 ; do
    echo "### Test Ex4-6-g.$i.$j"
    ../bin/carat/Conj_bravais -i "Ex4_g.$i" "Ex4_L.$j" > "g.$i.$j"
    echo "### Ex4-6-g.$i.$j return code $?"
  done
done
cat g.?.?

for i in 2 3 4 5 6 ; do
  for j in 2 3 4 ; do
    echo "### Test Ex4-7-pg.$i.$j"
    ../bin/carat/Extract -p "g.$i.$j" > "Ex4_pg.$i.$j"
    echo "### Ex4-7-pg.$i.$j return code $?"
    echo "### Test Ex4-7-cg.$i.$j"
    ../bin/carat/Extract -c "g.$i.$j" > "Ex4_cg.$i.$j"
    echo "### Ex4-7-cg.$i.$j return code $?"
  done
done
cat Ex4_pg.?.?
cat Ex4_cg.?.?

rm -f g.?.? Ex4_pg.?.? Ex4_cg.?.?
