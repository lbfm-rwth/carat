#include<typedef.h>
#include<base.h>
#include<graph.h>
#include<zass.h>
#include<tsubgroups.h>
#include<matrix.h>



/* ------------------------------------------------------------------ */
/* Berechne fast die Punktgruppe des affinen Normalisators von R:     */
/* nur Stab_{N_Gl_n(Z)(G)} (coz) ohne Punktgruppe von R               */
/* ------------------------------------------------------------------ */ 
/* R: Raumgruppe in Standard_affine_form ohne Translationen           */
/* P: Punktgruppe von R, P->gen und P->normal muessen                 */
/*    zusammen den Normalisator von P in Gl_n(Z) erzeugen.            */
/*    (P->zentr wird nicht beruecksichtigt)                           */
/* pres: Praesentation von R                                          */
/* anz: speichere die Anzahl der Matrizen dort                        */
/* ------------------------------------------------------------------ */
matrix_TYP **PoaN(bravais_TYP *R,
                  bravais_TYP *P,
		  matrix_TYP *pres,
		  int *anz)
{
   matrix_TYP **N, *coz, **H_G_Z;

   word *relator;

   int i;

   coz_TYP coz_info;


   /* Vorbereitungen */
   /* ============== */
   relator = (word *)calloc(pres->rows, sizeof(word));
   for (i = 0; i < pres->rows; i++){
      matrix_2_word(pres, relator + i, i);
   }
   coz = extract_c(R);

   /* calculate H^1(P, Q^n/Z^n) */
   H_G_Z = calculate_H1(P, relator, pres->rows);

   /* identify the cocyle */
   coz_info = identify_coz(P, coz, H_G_Z);

   /* get stabilizer */
   N = coz_info.Stab;
   coz_info.Stab = NULL;
   anz[0] = coz_info.Stab_no;

   /* Speicherfreigabe */
   /* ================ */
   free_mat(coz);
   for (i = 0; i < 3; i++)
      free_mat(H_G_Z[i]);
   free(H_G_Z);
   for (i = 0; i < pres->rows; i++)
      wordfree(relator + i);
   free(relator);
   free_coz_TYP(coz_info);

   return(N);
}
