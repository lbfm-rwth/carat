#include "typedef.h"
#include "name.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include "bravais.h"
#include "contrib.h"
#include "sort.h"
#include "datei.h"

/*************************************************************************
@ 
@ matrix_TYP *q_class_inf (bravais_TYP *G,
@                          database *database,
@                          char *qclass_name,
@                          char *symbol,
@                          bravais_TYP **OUT,
@                          matrix_TYP **PRES,
@                          int transformation)
@
@ SIDEEFECTS: G->gen[i] might be checked via Check_mat.
@             G->order and G->divisors WILL be set to the correct value.
**************************************************************************/
matrix_TYP *q_class_inf (bravais_TYP *G,
                         database *database,
                         char *qclass_name,
                         char *symbol,
                         bravais_TYP **OUT,
                         matrix_TYP **PRES,
                         int transformation)
{

  bravais_TYP *H;

  matrix_TYP  *DATABASE_MAT,
              *MATG,
              *ERG = NULL;  

  char filetmp[1024], format[1024];

  int i,
      orderG = 0,
      possible[1024],
      no_possible = 0;

  MATG = compute_q_matrix (G);

  for (i=0;i<MATG->cols;i++)
     orderG += MATG->array.SZ[0][i];

  G->order = orderG;
  factorize_new(G->order,G->divisors);

  get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/char.%s");
  for (i=0;i<database->nr;i++){
      if (orderG == database->entry[i].order &&
          G->dim == database->entry[i].degree && 
          MATG->rows == database->entry[i].no_idem + 9 &&
          MATG->cols == database->entry[i].no_conclass ){

          sprintf(filetmp, format,
                          G->dim,database->entry[i].symbol,
                          orderG,database->entry[i].discriminant,
                          database->entry[i].abbreviation);

          DATABASE_MAT = get_mat(filetmp);
          if (mat_comp(MATG,DATABASE_MAT) == 0){
               possible[no_possible] = i;
               no_possible++;
          }

          free_mat(DATABASE_MAT);
      }
  }


  /* there are 5 pairs of groups where the characteristic matrix does
     not decide the Q-equivalence. In these cases, get the hands dirty */
  get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/%s");
  while (no_possible > 1) {
     i = possible[no_possible-1];
     sprintf(filetmp, format,
                     G->dim,database->entry[i].symbol,
                     orderG,database->entry[i].discriminant,
                     database->entry[i].abbreviation);

     H = get_bravais(filetmp);
     ERG = suche_kand (G, H);

     if (ERG == 0){
        no_possible--;
        free_bravais(H);
     }
     else{
        possible[0] = i;
        no_possible = 1;
        if (OUT){
           *OUT = H;
        }
        else{
           free_bravais(H);
        }
     }
  }

  /* just for sanity */
  if (no_possible == 0){
    fprintf(stderr,"This group does not appear in the catalog of Q-classes.\n");
    fprintf(stderr,"Please verify that CARAT has been installed properly.\n");
    fprintf(stderr,"If so, please send a bug-report to\n");
    fprintf(stderr,"    carat@momo.math.rwth-aachen.de\n");
    fprintf(stderr,"including the forthcoming output.\n");
    put_bravais(G,NULL,NULL);
    exit(4);
  }

  get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/%s"); 
  if (OUT || (transformation && !ERG) ){
     i = possible[0];
     sprintf(filetmp, format,
                     G->dim,database->entry[i].symbol,
                     orderG,database->entry[i].discriminant,
                     database->entry[i].abbreviation);

     H = get_bravais(filetmp);

     if (OUT){
        *OUT = H;
     }

  }

  get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/pres.%s"); 
  if (PRES){
    i = possible[0];
     sprintf(filetmp, format,
                     G->dim,database->entry[i].symbol,
                     orderG,database->entry[i].discriminant,
                     database->entry[i].abbreviation);

     *PRES = get_mat(filetmp);
  }

  /* now we have a group which fits to G, lets see what we want to
     return as additional information */
  if (transformation && !ERG){

     ERG = suche_kand (G, H);

     if (OUT){
        *OUT = H;
     }
     else{
        free_bravais(H);
     }

     if (!ERG){
       fprintf(stderr,"This group does not appear in the catalog of Q-classes.\n");
       fprintf(stderr,"Please verify that CARAT has been installed properly.\n");
       fprintf(stderr,"If so, please send a bug-report to\n");
       fprintf(stderr,"    carat@momo.math.rwth-aachen.de\n");
       fprintf(stderr,"including the forthcoming output.\n");
       put_bravais(G,NULL,NULL);
       exit(4);
     }
  }

  i = possible[no_possible-1];
  if (qclass_name)
     sprintf(qclass_name,"%s",database->entry[i].abbreviation);

  if (symbol)
     sprintf(symbol,"%s",database->entry[i].symbol);

  free_mat(MATG);

  return ERG;

}

