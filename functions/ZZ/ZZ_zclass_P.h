extern int deal_with_ZCLASS _ZZ_P_PROTO_(( ZZ_data_t *data,
                       ZZ_tree_t *tree,
                       ZZ_node_t *father,
                       ZZ_node_t *new));
extern int orbit_under_normalizer _ZZ_P_PROTO_((
     ZZ_data_t *data,
     ZZ_tree_t *tree,
     ZZ_node_t *father,
     ZZ_node_t *new,
     int ii, 
     int jj,
     QtoZ_TYP *inzidenz,
     int *nr,
     ZZ_node_t **nnn));
extern matrix_TYP *special_deal_with_zclass _ZZ_P_PROTO_((ZZ_tree_t *tree,
                                                          bravais_TYP *group,
                                                          int *nr));
extern int in_bahn _ZZ_P_PROTO_((matrix_TYP *lattice,
                                 ZZ_node_t *father,
                                 int *i));

extern matrix_TYP *konjugierende _ZZ_P_PROTO_((matrix_TYP *L,
                                               bravais_TYP *G,
                                               ZZ_node_t *n));
