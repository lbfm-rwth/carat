#include"typedef.h"
#include"bravais.h"
#include"base.h"
#include"voronoi.h"
#include"matrix.h"
#include"getput.h"
#include"symm.h"
#include"autgrp.h"
#include"reduction.h"
#include"datei.h"



/* -----------------------------------------------------------
Berechne Normalisator von H:
H: Gruppe                                (WIRD GEAENDERT!!!!!!!!!!!!!!)
Gtr: transponierte von H, kann NULL sein
A: G-perfekte Form, kann NULL sein
prime: Primzahl, kann 0 sein
b_option: Normalisator der Bravaisgruppe von H wird berechnet
o_option: Gib Ergebnis auf Bildschirm aus
------------------------------------------------------------- */
void normalisator(bravais_TYP *H,
                  bravais_TYP *Gtr,
						matrix_TYP *A,
						int prime,
                  boolean b_option,
						boolean o_option)
{
   bravais_TYP *G;
	bravais_TYP *GB; /* the bravais group of G */

   matrix_TYP *tmp,
              *tmp2;

   voronoi_TYP **V;

   int Vanz,
       i;
       
   boolean gtr_flag = FALSE,
           a_flag = FALSE;



   /* throw away the normalizer of H */
   for (i = 0; i < H->normal_no; i++){
      free_mat(H->normal[i]);
   }
   if (H->normal != NULL && H->normal_no > 0)
      free(H->normal);
   H->normal = NULL;
   H->normal_no = 0;

   if (b_option){
      G = copy_bravais(H);
   }
   else{
      G = bravais_group(H, FALSE);
   }

   /* let's see whether we already got the formspace */
   if (G->form == NULL){
      G->form = formspace(G->gen, G->gen_no, 1, &G->form_no);
   }

   /* transposed group */
   if (Gtr == NULL){
      Gtr = tr_bravais(G,1,FALSE);
      gtr_flag = TRUE;
   }

   /* G-perfect form */
   if (A == NULL){
      /* firstly calculate an positive definite G-invariant form */
      tmp2 = init_mat(G->dim, G->dim, "1");
      tmp = rform(G->gen, G->gen_no, tmp2, 101);
      free_mat(tmp2);

      /* now calculate the trace bifo */
      tmp2 = trace_bifo(G->form ,Gtr->form, G->form_no);
      A = first_perfect(tmp, G, Gtr->form, tmp2, &Vanz);
      free_mat(tmp2);
      free_mat(tmp);
      
      a_flag = TRUE;
   }
   if (prime == 0)
      prime = 1949;

	GB = bravais_group(G,TRUE);
   V = normalizer(A, GB, Gtr, prime, &Vanz);
	G->normal = GB->normal; G->normal_no = GB->normal_no;
	GB->normal = NULL; GB->normal_no = 0;
	free_bravais(GB);

   /* now we got G and it's normalizer, so see what we can do with H */
   if (b_option){
      /* very easy in this case */
      H->normal = G->normal;
      H->normal_no = G->normal_no;
      G->normal = NULL;
      G->normal_no = 0;
   }
   else{
      H->normal = normalizer_in_N(H, G, &H->normal_no, FALSE);
   }

   if(o_option){
      put_bravais(H, NULL, "group with complete normalizer");
      put_bravais(Gtr, NULL, "Transposed group");
      for(i = 0; i < Vanz; i++)
         put_voronoi(V[i]);
   }

   /* cleaning up the memory */
   for(i=0;i<Vanz;i++){
      clear_voronoi(V[i]);
      free(V[i]);
   }
   free(V);
   free_bravais(G);
   if (a_flag)
      free_mat(A);
   if (gtr_flag)
      free_bravais(Gtr);
}

