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

#define DATABASE_NAME TOPDIR "/tables/qcatalog/data"

int SFLAG;
int INFO_LEVEL;
boolean GRAPH = FALSE;



bravais_TYP *get_qclass_by_name(char *name,
                                matrix_TYP **PRES,
                                int dim)
{

   int found = FALSE,
       i;

   database *database;

   char filename[1024];

   bravais_TYP *G;

   while (!found){
      database = load_database (DATABASE_NAME,dim);

      i = 0;
      while (!found && i<database->nr ){
         found = (strcmp(database->entry[i].abbreviation,name) == 0);
         i++;
      }

      if (found){
         i--;

         sprintf(filename,"%s/tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/%s",
                          TOPDIR,dim,database->entry[i].symbol,
                          database->entry[i].order,
                          database->entry[i].discriminant,name);
         G = get_bravais(filename);

         sprintf(filename,"%s/tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/pres.%s",
                          TOPDIR,dim,database->entry[i].symbol,
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


int main (int argc, char *argv[])
{
  
  bravais_TYP *R,
              *DATAQ,
              *DATAZ;

  matrix_TYP *PRES;

  MP_INT aff_name;

  char qname[1024],
      *affstring;


  int zname[2],
      i;

  char comment[1024];

  read_header (argc, argv);

    
  if (is_option('h'))
    INFO_LEVEL = optionnumber('h');
  
  if (INFO_LEVEL == 8)
    SFLAG = 1; 
  
  if ((is_option('h') && INFO_LEVEL != 8) || FILEANZ > 0)
    {
      printf("Usage: %s [-h] [-c] [-i]\n",argv[0]);
      printf("\n");
      printf("The program enables you to construct a space group corresponding\n");
      printf("to a name as given by CARAT. The program will read the various\n");
      printf("components of the name for the desired space group from stdin.\n");
      printf("It will then output generators for the space group, which are\n");
      printf("understood to generate the space group together with Z^n.\n");
      printf("Please note: not every name is a valid name!\n");
      printf("\n");
      printf("Options:\n");
      printf("-h     : gives you this help.\n");
      printf("-c     : do not check the name given, ie. verify that it is a valid name.\n");
      printf("         WARNING: this could lead to a wrong name in the header of the\n");
      printf("         resulting group\n");
      printf("-i     : ignore that the name is invalid, and give a space group at least\n");
      printf("         in the desired Q-class, and if exists one in the desired Z-class.\n");
      printf("         If given without -c the resulting group will also indicate the valid\n");
      printf("         name.\n");
      printf("Note: The Q-classes in               correspond to the Q-classes in degree\n");
      printf("   min.1 max.1                                                     1\n");
      printf("   min.2-5,group.1-4, max.2-3                                      2\n");
      printf("   min.6-14,group.5-25,max.4-5                                     3\n");
      printf("   min.15-57,group.26-205,max.6-9                                  4\n");
      printf("   min.58-169,group.206-1042,max.10-15                             5\n");
      printf("   min.170-667,group.1043-7636,max.16-27                           6\n");
      printf("\n");
      printf("Cf. Name, Q_catalog, QtoZ, Extensions.\n");

      if (FILEANZ == 0)
         exit(0);
      else
         exit(31);
    }


  fprintf(stderr,"qname:\n");
  scanf("%s",qname);

  fprintf(stderr,"zname: \n");
  scanf("%d %d",zname,zname+1);

  fprintf(stderr,"affname: \n");
  mpz_init(&aff_name);
  mpz_inp_str(&aff_name,stdin,10);


  if (is_option('c')){
     if (is_option('i')){
        i = 2;
     }
     else{
        i = 0;
     }
  }
  else{
     if (is_option('i')){
        i = 3;
     }
     else{
        i = 1;
     }
  }

  DATAQ = get_qclass_by_name(qname,&PRES,1);
  DATAZ = get_zclass_by_name(DATAQ,zname,zname+1,is_option('i'));
  R = get_affine_class_by_name(DATAZ,PRES,&aff_name,i);


  affstring = mpz_get_str (NULL,10,&aff_name);    
  sprintf(comment,"standart group with name %s %d %d %s",
                   qname,zname[0],zname[1],affstring);
  put_bravais(R,NULL,comment);
  free(affstring);
  free_mat(PRES);
  free_bravais(DATAZ);
  free_bravais(DATAQ);
  free_bravais(R);
  cleanup_prime();
  mpz_clear(&aff_name);

  if (INFO_LEVEL == 8) pointer_statistics(0,0);
  
  exit(0);
}
