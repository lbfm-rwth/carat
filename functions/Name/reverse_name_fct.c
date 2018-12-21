#include "ZZ.h"
#include "typedef.h"
#include "getput.h"
#include "name.h"
#include "bravais.h"
#include "datei.h"
#include "matrix.h"
#include "voronoi.h"
#include "autgrp.h"
#include "symm.h"
#include "base.h"
#include "zass.h"
#include "gmp.h"
#include "longtools.h"


bravais_TYP *get_qclass_by_name(char *name,
                                matrix_TYP **PRES,
                                int dim)
{

   int found = FALSE,
       i;

   database *database;

   char filename[1024], dbname[1024], format1[1024], format2[1024];

   bravais_TYP *G;

   get_data_dir(dbname,  "tables/qcatalog/data");
   get_data_dir(format1, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/%s");
   get_data_dir(format2, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/pres.%s");
   while (!found){
      database = load_database (dbname, dim);

      i = 0;
      while (!found && i<database->nr ){
         found = (strcmp(database->entry[i].abbreviation,name) == 0);
         i++;
      }

      if (found){
         i--;

         sprintf(filename, format1,
                          dim,database->entry[i].symbol,
                          database->entry[i].order,
                          database->entry[i].discriminant,name);
         G = get_bravais(filename);

         sprintf(filename, format2,
                          dim,database->entry[i].symbol,
                          database->entry[i].order,
                          database->entry[i].discriminant,name);
         *PRES = get_mat(filename);
      }
      else if (dim == 6){
         fprintf(stderr,"group with this name not in the catalog of Q-classes\n");
         fprintf(stderr,"name was %s\n",name);
         exit(4);
      }

      free_database (database);
      dim++;
   }

   return G;

}

bravais_TYP *get_zclass_by_name(bravais_TYP *G,
                                int *first,
                                int *second,
                                int ignore)
{

   bravais_TYP **QCLASS,
                *H,
                *RES,
               **TMP;

   int number,
       i;


   /* calculate the space of invariant forms of G, we'll need it */
   if (G->form_no == 0)
      G->form = formspace(G->gen,G->gen_no,1,&G->form_no);

   QCLASS = q2z(G,&number,TRUE,NULL,TRUE);


   if (number < *first){
      if (ignore){
         *first = number;
      }
      else{
         fprintf(stderr,"there is no zclass with this number this Q-class\n");
         exit(4);
      }
   }

   for (i=0;i<number;i++){
      if (i+1 != *first){
         free_bravais(QCLASS[i]);
         QCLASS[i] = NULL;
      } 
   }

   H = QCLASS[*first - 1];
   free(QCLASS);

   if (H->zentr_no < *second){
      if (ignore){
         *second = H->zentr_no;
      }
      else{
         fprintf(stderr,"there is no zclass with this number this Q-class\n");
         exit(4);
      }
   }

   for (i=0;i<H->zentr_no;i++){
      if (i+1 != *second){
         free_mat(H->zentr[i]);
         H->zentr[i] = NULL;
      }
   }

   TMP = get_groups(&H,1,&i);
   RES = TMP[*second - 1];

   free_bravais(H);
   free(TMP);

   if (RES == NULL){
      fprintf(stderr,"error in get_zclass_by_name\n");
      exit(4);
   }

   return RES;

}


bravais_TYP *split_extension(bravais_TYP *G)
{

   int i;

   bravais_TYP *R;

   R = init_bravais(G->dim+1);
   R->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
   for (i=0;i<G->gen_no;i++){
      R->gen[i] = copy_mat(G->gen[i]);
      extend(R->gen[i]);
   }
   R->gen_no = G->gen_no;

   return R;

}

bravais_TYP *get_affine_class_by_name(bravais_TYP *G,
                                     matrix_TYP *PRES,
                                     MP_INT *aff_name,
                                     int check)
{

   bravais_TYP *R;

   matrix_TYP **COH,
               *coz_mat,
              **matinv,
              **Y;

   word *relator;

   MP_INT standart_name,
          coho_size;

   int i,
       j,
       den;

   long dim;


   if (mpz_cmp_si(aff_name,0) == 0){
     /* handle the split extension case seperately */
     return split_extension(G);
   }

   /* do the cohomology calculations */
   relator = (word *) calloc(PRES->rows,sizeof(word));
   for (i=0;i<PRES->rows;i++){
     matrix_2_word(PRES,relator+i,i);
   }


   matinv = (matrix_TYP **) calloc(G->gen_no , sizeof(matrix_TYP *));
   COH = cohomology(&dim,G->gen,matinv,relator,G->gen_no,PRES->rows);

   if (dim > 0){
      fprintf(stderr,"error in get_affine_class_by_name\n");
      exit(4);
   }

   if (COH[0]->cols > 0){

      mpz_init(&coho_size);
      mpz_set_si(&coho_size,1);
      for (i=0;i<COH[1]->rows;i++)
         if (COH[1]->array.SZ[i][i]){
            mpz_mul_ui(&coho_size,&coho_size,
                      (unsigned long)COH[1]->array.SZ[i][i]);
         }

      mpz_init(&standart_name);
      mpz_mod(&standart_name,aff_name,&coho_size);

      coz_mat = reverse_valuation(&standart_name,COH[1]);

      R = space_group_from_matrix(G,coz_mat,COH[0],COH[1]);

      if (check == 1 || check == 3){
	 /* test whether the name we have got so far is the
           name CARAT would give */
         den = 1;
         for (i=0;i<G->gen_no;i++)
             den = KGV(R->gen[i]->kgv,den);
         free_mat(coz_mat);
         coz_mat = init_mat(G->dim*G->gen_no,1,"i");
         for (i=0;i<G->gen_no;i++)
	   for (j=0;j<G->dim;j++)
              coz_mat->array.SZ[j + G->dim * i][0] = den / R->gen[i]->kgv *
                    R->gen[i]->array.SZ[j][G->dim];
         coz_mat->kgv = den;

         Y = identify(COH[0],COH[1],COH[2],G,&coz_mat,&standart_name,1,
                     TRUE,NULL,NULL);

         /* now standart_name will contain the name CARAT would assign
            to this group. check ist now */
         if (mpz_cmp(&standart_name,aff_name) != 0){
	    if (check == 3){
               mpz_set(aff_name,&standart_name);
            }
            else{
               fprintf(stderr,"the given name is not CARAT`s name\n");
               exit(4);
            }
         }

         free_mat(Y[0]); free(Y);
      }

      free_mat(coz_mat);
      mpz_clear(&standart_name);
      mpz_clear(&coho_size);
   }
   else{
     /* things are easy here, there is only the split extension */
     R = split_extension(G);

     if (check == 1 && mpz_cmp_si(aff_name,0) != 0){
        fprintf(stderr,"wrong name for the space group\n");
        exit(4);
     }
     else{
        mpz_set_si(aff_name,0);
     }

   }

   for (i=0;i<G->gen_no;i++){
       if (matinv[i] != NULL)
          free_mat(matinv[i]);
   }
   free(matinv);
   for (i=0;i<3;i++) free_mat(COH[i]); free(COH);
   for (i=0;i<PRES->rows;i++) wordfree(relator+i);
   free(relator);


   return R;
}


bravais_TYP *reverse_name(char qname[1024],
                          int zname[2],
			  MP_INT aff_name,
			  int i,
			  boolean iflag,
			  char **affstring)
{
  bravais_TYP *R,
              *DATAQ,
              *DATAZ;

  matrix_TYP *PRES;



  DATAQ = get_qclass_by_name(qname, &PRES, 1);
  DATAZ = get_zclass_by_name(DATAQ, zname, zname+1, iflag);
  R = get_affine_class_by_name(DATAZ, PRES, &aff_name, i);

  affstring[0] = (char * ) malloc( mpz_sizeinbase( &aff_name , 10) + 2 );
  mpz_get_str (affstring[0] , 10, &aff_name);


  free_mat(PRES);
  free_bravais(DATAZ);
  free_bravais(DATAQ);

  cleanup_prime();
  mpz_clear(&aff_name);

  return(R);
}



