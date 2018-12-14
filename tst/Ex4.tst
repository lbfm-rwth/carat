
# compared to Ex4_gn, the result has a different set of normalizer generators,
# but the normalizer is the same
../bin/carat/Normalizer Ex4_g

../bin/carat/Vector_systems Ex4_gn

../bin/carat/Extract -r Ex4_gn Ex4_out

../bin/carat/Tr_bravais Ex4_gn

../bin/carat/Sublattices -b Ex4_gn_tr

for i in 2 3 4 5 6 ; do
  for j in 2 3 4 ; do
    ../bin/carat/Conj_bravais -i "Ex4_g.$i" "Ex4_L.$j" > "g.$i.$j"
  done
done
cat g.?.?

for i in 2 3 4 5 6 ; do
  for j in 2 3 4 ; do
    ../bin/carat/Extract -p "g.$i.$j" > "Ex4_pg.$i.$j"
    ../bin/carat/Extract -c "g.$i.$j" > "Ex4_cg.$i.$j"
  done
done
cat Ex4_pg.?.?
cat Ex4_cg.?.?

rm -f g.?.? Ex4_pg.?.? Ex4_cg.?.?
