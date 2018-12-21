/* Dieser Quellcode wurde name.c entnommen! */

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
#include "tsubgroups.h"


/* ------------------------------------------------------------------- */
/* Berechne den Namen einer Raumgruppe R in CARAT-Notation und gebe    */
/* diesen aus. P sei die Punktgruppe von R mit korrektem Formenraum.   */
/* ------------------------------------------------------------------- */
CARATname_TYP name_fct(bravais_TYP *R,
                       database *database)
{

  bravais_TYP *P,
              *Rnew,
              *RC,
              *DATAQ = (bravais_TYP *) 1,
              *DATAZ;

  matrix_TYP *T,
             *TI,
             *TZ,
             *PRES;

  MP_INT aff_name;

  char symb[1024];

  CARATname_TYP Name;



  P = point_group(R, 2);
  T = q_class_inf (P,database,Name.qname,symb,&DATAQ,&PRES,FALSE);
  TZ = z_class_inf(P,DATAQ,&DATAZ,Name.zname);
  if (DATAQ->order == 0)
     Name.order = 1;
  else
     Name.order = DATAQ->order;
  free_bravais(DATAQ);


  extend(TZ);
  Rnew = konj_bravais(R,TZ);
  T = mat_inv(TZ); free_mat(TZ); TZ = T; T = NULL;
  mpz_init(&Name.aff_name);

  /*
  if (is_option('o'))
     T = aff_class_inf(Rnew,DATAZ,PRES,&aff_name,&RC);
  else */
     T = aff_class_inf(Rnew,DATAZ,PRES,&Name.aff_name,NULL);

  Check_mat(T);
  Check_mat(TZ);
  TI = long_mat_inv(T);
  mat_muleq(TZ,TI);
  Name.trafo = TZ;

  /*
  if (is_option('o')){
     sprintf(comment,"standard group for %s",FILENAMES[0]);
     put_bravais(RC,NULL,comment);
     free_bravais(RC);
  }

  if (is_option('c')){	// Oliver: 10.04.2002
     printf("%s-%d.%d-", Name.qname, Name.zname[0], Name.zname[1]);
     mpz_out_str(stdout, 10, &aff_name);
     printf("\n");
  }
  else{
     printf("qname: %s ",Name.qname);
     printf("zname: %d %d ",Name.zname[0],Name.zname[1]);
     printf("aff_name: "); mpz_out_str(stdout,10,&aff_name); printf("\n");
  }
  */

  free_mat(PRES);
  free_mat(T);
  free_mat(TI);
  free_bravais(Rnew);
  free_bravais(P);
  free_bravais(DATAZ);
  cleanup_prime();

  return(Name);
}
