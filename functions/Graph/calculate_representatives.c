/* last change: 24.05.2002 by Oliver Heidbuechel */


#include <ZZ.h>
#include <typedef.h>
#include <presentation.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <datei.h>
#include <graph.h>
#include <gmp.h>
#include <zass.h>
#include <tools.h>
#include <longtools.h>
#include <sort.h>


/* ================================================================= */
/* Der Quellcode wurde normalop.c (identify) entnommen!!!            */
/* Stelle G->normal auf H^1(G, Q^n/Z^n) dar.                         */
/* ================================================================= */
/* G: Gruppe							     */
/* H_G_Z: H^1(G, Q^n/Z^n)                                            */
/* ================================================================= */
static matrix_TYP **stelle_normalisator_dar(bravais_TYP *G,
                                            matrix_TYP **H_G_Z)
{
  matrix_TYP **N;

  int i;

  /* put_bravais(G, 0, 0); */

  N = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));
  for (i = 0; i < G->normal_no; i++){
     N[i] = normalop(H_G_Z[0], H_G_Z[1], H_G_Z[2], G, G->normal[i], i == (G->normal_no - 1));
     if (INFO_LEVEL & 16){
        put_mat(N[i], NULL, "N[i]", 2);
     }
  }

  return(N);
}

/*
Uebersetzungstabelle:

cozykle:	H_G_Z[0]
D:		H_G_Z[1]
R:		H_G_Z[2]

*/



/* ================================================================= */
/* Der Quellcode wurde dem Programm Extensions entnommen!!!          */
/* Berechne die Ordnung von H^1(G,Q^n/Z^n)!                          */
/* ================================================================= */
/* H_G_Z: H^1(G, Q^n/Z^n)                                            */
/* ================================================================= */
static int CohomSize(matrix_TYP **H_G_Z)
{
  int n;

  MP_INT cohom_size;


  if (H_G_Z[0]->cols <1){
     return(1);
  }
  else {
     cohom_size = cohomology_size(H_G_Z[1]);
     n = mpz_get_ui(&cohom_size);
     mpz_clear(&cohom_size);
     return(n);
  }
}




/* ================================================================= */
/* Berechne den Stabilisator von elem + Ker phi                      */
/* ================================================================= */
/* H1_mod_ker: Daten ueber H^1 / ker (phi)                           */
/* S1: Stabilisator                                                  */
/* new_rep: representation from S1 on H^1(G,Q^n/lattice) / Ker(phi)  */
/* S1_word_no: Anzahl der Elemente von new_rep                       */
/* elem: Repraesentant von elem + Ker phi                            */
/* H_GL_Z: H^1(GL, Q^n/Z^n) isom. zu H^1(G, Q^n/L)                   */
/* counter: hier wird die Anzahl der Stabilisatorerzeuger gespeichert*/
/* ================================================================= */
static matrix_TYP **stab_preimage_Ker(H1_mod_ker_TYP H1_mod_ker,
                                      matrix_TYP **S1,
                                      matrix_TYP **new_rep,
				      int S1_word_no,
				      matrix_TYP *elem,
				      matrix_TYP **H_GL_Z,
				      int *counter)
{
   int reps_no, *length, m, k,
       ***WORDS, *WORDS_no;

   matrix_TYP ***reps, **S_ksi, *temp, **S1_inv, *std;



   /* Berechne Standardform von elem als Element von H^1/Ker */
   /*if (elem == NULL){
      temp = init_mat(H1_mod_ker.i->rows, 1, "");
   }
   else{*/
      temp = mat_mul(H1_mod_ker.i, elem);
   /*}*/
   for (k = 0; k < temp->rows; k++){
      temp->array.SZ[k][0] %= H1_mod_ker.D->array.SZ[k][k];
      if (temp->array.SZ[k][0] < 0)
         temp->array.SZ[k][0] += H1_mod_ker.D->array.SZ[k][k];
   }
   std = init_mat(H1_mod_ker.erz_no, 1, "");
   for (k = 0; k < H1_mod_ker.erz_no; k++){
      std->array.SZ[k][0] = temp->array.SZ[k + H1_mod_ker.D_first][0];
   }
   free_mat(temp);

   /* Bahn von elem sowie dessen Stabilisator berechnen */
   reps = H1_mod_ker_orbit_alg(H1_mod_ker, new_rep, S1_word_no, &reps_no,
                              &length, &WORDS, &WORDS_no, std);
   free_mat(std);
   if (reps_no != 2 || length[1] < 1){
      fprintf(stderr, "ERROR 1 in stab_preimage_Ker!\n");
      exit(78);
   }

   /* Berechne Stabilisator */
   S1_inv = (matrix_TYP **)calloc(S1_word_no, sizeof(matrix_TYP *));
   counter[0] = 0;
   S_ksi = (matrix_TYP **)calloc(WORDS_no[1], sizeof(matrix_TYP *));
   for (m = 0; m < WORDS_no[1]; m++){
      /* simplify the words */
      normalize_word(WORDS[1][m]);
      if (WORDS[1][m][0] == 1 && WORDS[1][m][1] < 0)
         WORDS[1][m][1] *= (-1);
      if (word_already_there(WORDS[1], m) == 0){
         S_ksi[counter[0]] = graph_mapped_word(WORDS[1][m], S1, S1_inv, H_GL_Z[1]);
         if (mat_search(S_ksi[counter[0]], S_ksi, counter[0], mat_comp) != -1){
            free_mat(S_ksi[counter[0]]);
         }
         else{
            mat_quicksort(S_ksi, 0, counter[0], mat_comp);
            counter[0]++;
         }
      }
   }

   /* clean */
   /* ===== */
   for (m = 0; m < S1_word_no; m++){
      if (S1_inv[m] != NULL)
         free_mat(S1_inv[m]);
   }
   free(S1_inv);
   for (m = 0; m < WORDS_no[1]; m++){
      free(WORDS[1][m]);
   }
   for (m = 0; m < length[1]; m++){
      free_mat(reps[1][m]);
   }
   free(reps[1]);
   free(WORDS[1]);
   free(reps);
   free(WORDS);
   free(WORDS_no);
   free(length);

   return(S_ksi);
}




/* ================================================================= */
/* Der Quellcode wurde subgroupgraph.c entnommen!!!                  */
/* Berechne Vertreter der Bahnen der Untergruppen unter dem          */
/* Stab des Gitters geschnitten mit dem Stabilisator des Cocykels    */
/* als Element von H^1(G, Q^n/L)!!!                                  */
/* ================================================================= */
/* G: Punktruppe                                                     */
/* GL: G konjugiert mit Gitter lattice                               */
/* H_G_Z: H^1(G, Q^n/Z^n)                                            */
/* H_GL_Z: H^1(GL, Q^n/Z^n) isom. zu H^1(G, Q^n/L)                   */
/* lattice: Gitter                                                   */
/* phi: H^1(G, Q^n/L) isom. zu H^1(GL,Q^n,Z^n) -> H^1(G,Q^n/Z^n)     */
/* H_1_mod_ker: Daten ueber H^1 / ker (phi)                          */
/* preimage: Element aus H^1(GL, Q^n/Z^n), welches auf Cozykel zur   */
/*           Raumgruppe abgebildet wird.                             */
/* kernel_mat: Matrix, die den Kern von phi beschreibt               */
/* orbit_no: speichere die Anzahl der Vertreter hier                 */
/* orbit_length: speichere die Laengen der Bahnen hier               */
/* ================================================================= */
matrix_TYP **calculate_representatives(bravais_TYP *G,
                                       bravais_TYP *GL,
				       matrix_TYP **H_G_Z,
				       matrix_TYP **H_GL_Z,
				       matrix_TYP *lattice,
                                       matrix_TYP *phi,
                                       H1_mod_ker_TYP H1_mod_ker,
				       matrix_TYP *preimage,
				       matrix_TYP *kernel_mat,
				       int *orbit_no,
				       int **orbit_length)
{
   matrix_TYP **S1 = 0,
              **lNli,
	      **Ninv,
	      **new_rep,
	      *invlattice,
	      **N_H_GL_Z,
	      *id,
	      **S_ksi,
              **kernel_gen,
              **kernel_elements,
	      **orbit_rep,
	      **rep;

   int l,
       norm_no, S1_word_no, cohomSize, 
       kernel_order, S_ksi_no,
       *kernel_list;

   rational eins;



   eins.z = eins.n = 1;


   /* calculate S1 = Stab_N(L) */
   /* ======================== */
   if (phi->cols > 0){
      norm_no = GL->normal_no;
      invlattice = mat_inv(lattice);
      id = init_mat(G->dim, G->dim, "1");
      Ninv = (matrix_TYP **)calloc(norm_no, sizeof(matrix_TYP *));
      lNli = (matrix_TYP **)calloc(norm_no, sizeof(matrix_TYP *));
      for (l = 0; l < norm_no; l++){
         lNli[l] = mat_kon(lattice, GL->normal[l], invlattice);
      }
      free_mat(invlattice);
      N_H_GL_Z = stelle_normalisator_dar(GL, H_GL_Z);
      S1 = calculate_S1(id, lNli, norm_no, &S1_word_no, N_H_GL_Z, Ninv, H_GL_Z[1]);

      for (l = 0; l < norm_no; l++){
         free_mat(lNli[l]);
	 free_mat(N_H_GL_Z[l]);
	 if (Ninv[l] != NULL)
	    free_mat(Ninv[l]);
      }
      free_mat(id);
      free(lNli);
      free(Ninv);
      free(N_H_GL_Z);
   }
   else{
      /* trivial part */
      S1_word_no = 0;
   }

   cohomSize = CohomSize(H_GL_Z);

   /* Kern */
   kernel_gen = col_to_list(kernel_mat);
   kernel_elements = (matrix_TYP **)calloc(cohomSize, sizeof(matrix_TYP *));
   kernel_list = aufspannen(cohomSize, kernel_elements, kernel_gen,
                            kernel_mat->cols, H_GL_Z[1], &kernel_order);


   if (H1_mod_ker.flag == 0 && H1_mod_ker.erz_no > 0){
      if (preimage == NULL){
         rep = orbit_ker(kernel_elements, kernel_order, H_GL_Z[1], S1,
                         S1_word_no, cohomSize, kernel_list, orbit_no, orbit_length);
      }
      else{
         /* representation from S1 on H^1(G,Q^n/lattice) / Ker(phi) */
         new_rep = new_representation(S1, S1_word_no, H1_mod_ker, H_GL_Z[1]);

         /* affine representation of the kernel */
         kernel_elements_2_affine(kernel_elements, cohomSize);

         /* calculate Stab_S1 (preimage + Kern) */
         S_ksi = stab_preimage_Ker(H1_mod_ker, S1, new_rep, S1_word_no, preimage, H_GL_Z, &S_ksi_no);

         /* calculate representatives */
         orbit_rep = orbit_ksi_plus_ker(preimage, kernel_elements, kernel_order,
                                        H_GL_Z[1], S_ksi, S_ksi_no, cohomSize,
	   			        kernel_list, orbit_no, orbit_length);
         rep = (matrix_TYP **)calloc(orbit_no[0], sizeof(matrix_TYP *));
         for (l = 0; l < orbit_no[0]; l++){
            rep[l] = mat_add(orbit_rep[l], preimage, eins, eins);
	    free_mat(orbit_rep[l]);
	 }
	 free(orbit_rep);
         for (l = 0; l < S_ksi_no; l++)
            free_mat(S_ksi[l]);
         free(S_ksi);
         for (l = 0; l < S1_word_no; l++)
            free_mat(new_rep[l]);
         free(new_rep);
      }
   }
   else{
      switch (H1_mod_ker.flag){
         case 0:
         case 2:
            /* orbit on 0 + Ker(phi) */
            rep = orbit_ker(kernel_elements, kernel_order, H_GL_Z[1], S1,
                            S1_word_no, cohomSize, kernel_list, orbit_no, orbit_length);
            break;
         case 1:
         case 3:
            rep = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
	    rep[0] = init_mat(0, 0, "");
	    orbit_length[0] = (int *)calloc(1, sizeof(int));
	    orbit_length[0][0] = 1;
	    orbit_no[0] = 1;
            break;
	 default:
            fprintf(stderr, "ERROR 1 in cacluate_representative!\n");
            exit(55);
      }
   }


   /* clean */
   /* ===== */
   if (phi->cols > 0){
      for (l = 0; l < S1_word_no; l++){
         free_mat(S1[l]);
       }
       free(S1);
   }
   for (l = 0; l < cohomSize; l++){
      if (kernel_elements[l] != NULL)
         free_mat(kernel_elements[l]);
   }
   free(kernel_elements);
   free(kernel_list);
   for (l = 0; l < kernel_mat->cols; l++){
      free_mat(kernel_gen[l]);
   }
   free(kernel_gen);


   return(rep);
}




/*
Uebersetzungstabelle:
=====================

data->Z[j]		GL
data->Z[i]		G
data->X[j]		H_GL_Z
data->X[i]		H_G_Z
data->N[j]		N_H_GL_Z


*/

