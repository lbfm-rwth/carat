extern int deal_with_ZCLASS(ZZ_tree_t *tree, ZZ_node_t *newnode);
extern int orbit_under_normalizer(
     ZZ_tree_t *tree,
     ZZ_node_t *father,
     ZZ_node_t *newnode,
     QtoZ_TYP *inzidenz,
     int *nr,
     ZZ_node_t **nnn);
extern matrix_TYP *special_deal_with_zclass(ZZ_tree_t *tree,
                                                          bravais_TYP *group,
                                                          int *nr);
extern int in_bahn(matrix_TYP *lattice,
                                 ZZ_node_t *father,
                                 int *i);

extern matrix_TYP *konjugierende(matrix_TYP *L,
                                               bravais_TYP *G,
                                               ZZ_node_t *n);
