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
#include "contrib.h"
#include "base.h"
#include "zass.h"
#include "gmp.h"
#include "longtools.h"

boolean GRAPH = FALSE;
int SFLAG;
int INFO_LEVEL;
int main (int argc, char *argv[])



{
  
  bravais_TYP *P,
              *R,
              *Rnew,
              *RC,
              *DATAQ = (bravais_TYP *) 1,
              *DATAZ;

  matrix_TYP *T,
             *TI,
             *TZ,
             *PRES;

  database *database;

  MP_INT aff_name;

  char qname[1024],
       symb[1024];


  int zname[2];

  char comment[1024], dbname[1024];

  read_header (argc, argv);

    
  if (is_option('h'))
    INFO_LEVEL = optionnumber('h');
  
  if (INFO_LEVEL == 8)
    SFLAG = 1; 
  
  if ((is_option('h') && INFO_LEVEL != 8) || FILEANZ == 0)
    {
      printf("Usage: %s file [-T] [-Z] [-o] [-M] [-c]\n",argv[0]);
      printf("\n");
      printf("file: bravais_TYP containing the space group R or the finite unimodular\n");
      printf("      group G.\n");
      printf("\n");
      printf("The program is used to give crystallographic groups a name, i.e.\n");
      printf("compute a string which is depends only on the crystallographic\n");
      printf("class of the group R/G, and determines it uniquely.\n");
      printf("CAUTION: the program assumes the translation lattice to be Z^n.\n");
      printf("If called with -Z, the program assumes file to contain a finite unimodular\n");
      printf("group G, and will output a name for the arithmetic class of G.\n");
      printf("Otherwise, it will assume file to contain generators for the space group\n");
      printf("R, and will output a name for the affine class of R.\n");
      printf("\n");
      printf("Note that the first part of the name for a space group R is exactly the\n");
      printf("name for the arithmetic class of the point group G of R.\n");
      printf("\n");
      printf("Options:\n");
      printf("-h    : gives this help.\n");
      printf("-T    : output a matrix transforming the given group R/G into CARAT's\n");
      printf("        representative.\n");
      printf("-o    : output CARAT's representative for this affine/arithmetic class.\n");
      printf("-Z    : assume file to contain a finite unimodular group, more details\n");
      printf("        see above.\n");
      printf("-M    : give short Hermann-Mauguin symbols to describe a group\n");
      printf("        isomorphic to the given only (has an effect only if the\n");
      printf("        degree of R is two or three)\n");
      printf("-c    : gives the CARAT name as\n");	// Oliver: 10.04.2002
      printf("        \"Q-class-name\"-\"Z-class-name\"(-\"name for the affine class\")\n");
      printf("\n");
      printf("Cf.: Q_catalog, QtoZ, Extensions, Symbol, Standard_affine_form.\n");

      if (FILEANZ == 0)
         exit(0);
      else
         exit(31);
    }

  R = get_bravais(FILENAMES[0]);

  if (is_option('Z')) {
     P = copy_bravais(R);
  }
  else{
     P = point_group(R,2);
  }

  get_data_dir(dbname, "tables/qcatalog/data");
  database = load_database (dbname, P->dim);

  T = q_class_inf (P,database,qname,symb,&DATAQ,&PRES,FALSE);

  TZ = z_class_inf(P,DATAQ,&DATAZ,zname);
  free_bravais(DATAQ);


  if (is_option('Z')){
     /* we are finished here */
     if (is_option('c')){	// Oliver: 10.04.2002
        printf("%s-%d.%d\n",qname,zname[0],zname[1]);
     }
     else{
        printf("qname: %s ",qname);
        printf("zname: %d %d\n",zname[0],zname[1]);
     }

     if (is_option('T')){
        sprintf(comment,"transformation matrix for %s",FILENAMES[0]);
        put_mat(TZ,NULL,comment,2);
     }

     if (is_option('o')){
        sprintf(comment,"standard group for %s",FILENAMES[0]);
        put_bravais(DATAZ,NULL,comment);
     }

     free_bravais(P);
     free_bravais(R);
     free_mat(TZ);
     free_bravais(DATAZ);
     free_database (database);
     cleanup_prime();

     if (INFO_LEVEL == 8) pointer_statistics(0,0);

     exit(0);
  }

  extend(TZ);
  Rnew = konj_bravais(R,TZ);

  T = mat_inv(TZ); free_mat(TZ); TZ = T; T = NULL;

  if (INFO_LEVEL & 4)
     put_bravais(Rnew,NULL,NULL);

  if (INFO_LEVEL & 4)
     put_mat(TZ,NULL,"transformation matrix",2);


  mpz_init(&aff_name);

  if (is_option('o'))
     T = aff_class_inf(Rnew,DATAZ,PRES,&aff_name,&RC);
  else
     T = aff_class_inf(Rnew,DATAZ,PRES,&aff_name,NULL);

  /* put_mat(TZ,0,"TZ",0);
  put_mat(T,0,"T",0); */

  if (is_option('T')){
     Check_mat(T);
     Check_mat(TZ);
     TI = long_mat_inv(T);
     mat_muleq(TZ,TI);
     put_mat(TZ,NULL,"transformation matrix",2);
     free_mat(TI);
  }

  if (is_option('o')){
     sprintf(comment,"standard group for %s",FILENAMES[0]);
     put_bravais(RC,NULL,comment);
     free_bravais(RC);
  }

  if (is_option('c')){	// Oliver: 10.04.2002
     printf("%s-%d.%d-", qname, zname[0], zname[1]);
     mpz_out_str(stdout, 10, &aff_name);
     printf("\n");
  }
  else{
     printf("qname: %s ",qname);
     printf("zname: %d %d ",zname[0],zname[1]);
     printf("aff_name: "); mpz_out_str(stdout,10,&aff_name); printf("\n");
  }

  if (is_option('M') && (R->dim == 4 || R->dim == 3)){
     display_HM_symbol(qname,zname[0],zname[1],&aff_name);
  }

  free_mat(PRES);
  free_mat(TZ);
  free_mat(T);
  free_bravais(P);
  free_bravais(Rnew);
  free_bravais(R);
  free_bravais(DATAZ);
  free_database (database);
  mpz_clear(&aff_name);
  cleanup_prime();

  if (INFO_LEVEL == 8) pointer_statistics(0,0);

  exit(0);
}
