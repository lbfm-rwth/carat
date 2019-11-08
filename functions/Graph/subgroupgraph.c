/* last change: 16.03.01 by Oliver Heidbuechel */



#include <typedef.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <graph.h>
#include <zass.h>
#include <datei.h>
#include <longtools.h>
#include <sort.h>
#include <presentation.h>


typedef struct{
   int **x;
   int ***l;
   int ***o;
} SOL_TYP;



/* ----------------------------------------------------------------------------- */
static SOL_TYP *init_SOL_TYP(int gitter_no,
                             int aff_i,
                             int aff_j,
			     boolean oflag)
{
   int i, j;

   SOL_TYP *sol;



   sol = (SOL_TYP *)calloc(gitter_no, sizeof(SOL_TYP));

   for (i = 0; i < gitter_no; i++){
      sol[i].x = (int **)calloc(aff_i, sizeof(int *));
      sol[i].l = (int ***)calloc(aff_i, sizeof(int **));
      if (oflag)
         sol[i].o = (int ***)calloc(aff_i, sizeof(int **));
      for (j = 0; j < aff_i; j++){
         sol[i].x[j] = (int *)calloc(aff_j, sizeof(int));
         sol[i].l[j] = (int **)calloc(aff_j, sizeof(int *));
	 if (oflag)
	    sol[i].o[j] = (int **)calloc(aff_j, sizeof(int *));
      }
   }

   return(sol);
}



/* ----------------------------------------------------------------------------- */
static void free_SOL_TYP(SOL_TYP *sol,
                         int gitter_no,
                         int aff_i,
                         int aff_j)
{
   int i, j, k;

   for (i = 0; i < gitter_no; i++){
      for (j = 0; j < aff_i; j++){
         for (k = 0; k < aff_j; k++){
            if (sol[i].l[j][k] != NULL)
               free(sol[i].l[j][k]);
            if (sol[i].o != NULL && sol[i].o[j] != NULL && sol[i].o[j][k] != NULL)
               free(sol[i].o[j][k]);
         }
         free(sol[i].l[j]);
         free(sol[i].x[j]);
	 if (sol[i].o != NULL && sol[i].o[j] != NULL)
	    free(sol[i].o[j]);
      }
      free(sol[i].l);
      free(sol[i].x);
      if (sol[i].o != NULL)
         free(sol[i].o);
   }
   free(sol);
}



/* ----------------------------------------------------------------------------- */
static boolean gibt_gitter_beitrag(int *gitter_no,
                                    int laenge,
                                    int no)
{
   int i;

   for (i = 1 ; i <= laenge; i++){
      if (gitter_no[i] == no)
         return (TRUE);
   }
   return (FALSE);
}



/* ----------------------------------------------------------------------------- */
/* test if two groups are equal                                                  */
/* ----------------------------------------------------------------------------- */
static int test_fkt(bravais_TYP *t,
                    bravais_TYP *G)
{
   int i, j, k, counter, flagge, ret;

   matrix_TYP **Gorbit, *tmp;


   Gorbit = (matrix_TYP **)calloc(G->order, sizeof(matrix_TYP *));
   Gorbit[0] = init_mat(t->dim, t->dim, "1");
   counter = 1;
   i = 0;

   while (counter < G->order){
      for (j = 0; j < G->gen_no; j++){
         tmp = mat_mul(Gorbit[i], G->gen[j]);
         flagge = 1;
         for (k = 0; k < counter && flagge != 0; k++){
            flagge = cmp_mat(tmp, Gorbit[k]);
         }
         if (flagge != 0){
            Gorbit[counter] = tmp;
            counter++;
         }
         else{
            free_mat(tmp);
         }
      }
      i++;
   }

   ret = 1;
   for (i = 0; i < t->gen_no; i++){
      flagge = 0;
      for (j = 0; j < G->order; j++){
         if (cmp_mat(Gorbit[j], t->gen[i]) == 0){
            flagge = 1;
            break;
         }
      }
      if (flagge == 0){
         ret = 0;
         break;
      }
   }

   for (i = 0; i < G->order; i++){
      free_mat(Gorbit[i]);
   }
   free(Gorbit);

   return(ret);
}



/* ----------------------------------------------------------------------------- */
/* calculate graph of incidences for k-subgroups for a given geometric class     */
/* beeing not 1 or -1!                                                           */
/* ----------------------------------------------------------------------------- */
matrix_TYP *subgroupgraph(Q_data_TYP *data,
                          boolean oflag)
{
   int i, j, k, l, m, q, flagge, reps_no,
       image_gen_no, S1_word_no, norm_no, kernel_order, orbit_no,
       *kernel_list, *orbit_length, *length,
       ***WORDS, *WORDS_no,
       i__, j__, counter, translanz, *kernel_factor,
       lattice_orbits_no, *lattice_orbits, *lattice_orbit_length,
       **gitter_no, *smallest, i__first, reps_flag,
       *primes, *exponent, *konj_wort;

   QtoZ_entry_TYP entry;

   matrix_TYP *erg,
              *id,
              *lattice,
              *invlattice,
              *diag,
              *phi,
              *kernel_mat,
             **kernel_gen,
             **kernel_elements,
             **image,
             **lNli,
             **S1 = 0,
             **S1_inv = 0,
             **S_ksi,
            ***dataNinv,
             **orbit_rep,
             **new_rep = 0,
            ***reps = 0,
              *coz_j,
              *coz_i,
              *cocycle,
             **translationen,
             **TRASH;

   word *relator;

   bravais_TYP *G;

   H1_mod_ker_TYP H1_mod_ker;

   SOL_TYP *zwischenerg;

   rational eins;



   eins.z = eins.n = 1;

   dataNinv = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   for (i = 0; i < data->Z_no; i++){
      dataNinv[i] = (matrix_TYP **)calloc(data->Z[i]->normal_no, sizeof(matrix_TYP *));
   }
   id = init_mat(data->G->dim, data->G->dim, "1");
   erg = init_mat(data->all, data->all, "");

   /* convert presentation into relator format */
   relator = (word *) calloc(data->pres->rows,sizeof(word));
   for (i = 0; i < data->pres->rows; i++){
      matrix_2_word(data->pres, relator + i, i);
   }

   for (i = 0; i < data->Z_no; i++){

      translationen = transl_aff_normal(data->Z[i]->gen, data->Z[i]->gen_no, &translanz);

      for (j = 0; j < data->Z_no; j++){
         entry = data->INZ->entry[i][j];
         zwischenerg = init_SOL_TYP(entry.anz, data->aff_no[i], data->aff_no[j], oflag);
         gitter_no = (int **)calloc(data->aff_no[i], sizeof(int *));
         for (k = 0; k < data->aff_no[i]; k++){
            gitter_no[k] = (int *)calloc(entry.anz + 1, sizeof(int));
         }
         kernel_factor = (int *)calloc(entry.anz, sizeof(int));
         primes = (int *)calloc(entry.anz, sizeof(int));
         exponent = (int *)calloc(entry.anz, sizeof(int));
         for (k = 0; k < entry.anz; k++){
            lattice = entry.lattice[k];

            /* calculate index of the sublattice */
            for (l = 0; l < lattice->rows; l++){
               if (entry.lsf[k]->array.SZ[l][l] != 1){
                  primes[k] = entry.lsf[k]->array.SZ[l][l];
                  exponent[k]++;
               }
            }

            invlattice = mat_inv(lattice);
            diag = matrix_on_diagonal(lattice, data->G->gen_no);

            /* calculate information about Ker phi | B^1 */
            TRASH = kernel_factor_fct(translationen, translanz, data->Z[i]->gen_no, lattice,
                                      &kernel_factor[k]);
            for (l = 0; l < kernel_factor[k]; l++){
               free_mat(TRASH[l]);
            }
            free(TRASH);

            /* test if lattice^(-1) data->Z[i]->gen[l] lattice == data->Z[j]->gen[l] */
            flagge = 1;
            G = konj_bravais(data->Z[i], invlattice);
            if (GRAPH_DEBUG){
               if (!test_fkt(G, data->Z[j])){
                  fprintf(stderr, "ERROR 1 in subgroupgraph!\n");
                  exit(8);
               }
            }
            for (l = 0; l < G->gen_no; l++){
               if (cmp_mat(data->Z[j]->gen[l], G->gen[l]) != 0){
                  flagge = 0;
                  break;
               }
            }

            /* calculate phi: H^1(Z[i], Q^n/lattice) -> H^1(Z[i], Q^n/Z^n) */
            if (flagge == 0){
               cocycle = H1_of_standard_to_GL(G, data->Z[j], data->X[j]);
               calculate_phi(diag, cocycle, data->X[i], data->X[j], data->X_2_inv[i], &phi,
                             &kernel_mat, &image, &image_gen_no, &H1_mod_ker);
               free_mat(cocycle);
            }
            else{
               calculate_phi(diag, data->X[j][0], data->X[i], data->X[j], data->X_2_inv[i], &phi,
                             &kernel_mat, &image, &image_gen_no, &H1_mod_ker);
            }
            kernel_gen = col_to_list(kernel_mat);
            kernel_elements = (matrix_TYP **)calloc(data->coho_size[j], sizeof(matrix_TYP *));
            kernel_list = aufspannen(data->coho_size[j], kernel_elements, kernel_gen,
                                        kernel_mat->cols, data->X[j][1], &kernel_order);
            free_bravais(G);

            /* calculate S1 = Stab_N(L)*/
            if (phi->cols > 0){
               norm_no = data->Z[j]->normal_no;
               lNli = (matrix_TYP **)calloc(norm_no, sizeof(matrix_TYP *));
               for (l = 0; l < norm_no; l++){
                  lNli[l] = mat_kon(lattice, data->Z[j]->normal[l], invlattice);
               }

               S1 = calculate_S1(id, lNli, norm_no, &S1_word_no,
                                 data->N[j], dataNinv[j], data->X[j][1]);

               S1_inv = (matrix_TYP **)calloc(S1_word_no, sizeof(matrix_TYP *));
               for (l = 0; l < norm_no; l++){
                  free_mat(lNli[l]);
               }
               free(lNli);
            }
            else{
               /* trivial part */
               S1_word_no = 0;
            }

            if (H1_mod_ker.flag == 0 && H1_mod_ker.erz_no > 0){
               /* representation from S1 on H^1(G,Q^n/lattice) / Ker(phi) */
               new_rep = new_representation(S1, S1_word_no, H1_mod_ker, data->X[j][1]);

               /* orbit on H^1(G,Q^n/lattice) / Ker(phi) */
               reps = H1_mod_ker_orbit_alg(H1_mod_ker, new_rep, S1_word_no, &reps_no, &length,
                                           &WORDS, &WORDS_no, NULL);

               /* l = 0 (d.h. Untergruppen der zerfallenden Gruppe) */
	       /* orbit on 0 + Ker(phi): l = 0 */
               orbit_rep = orbit_ker(kernel_elements, kernel_order, data->X[j][1], S1,
                                     S1_word_no, data->coho_size[j], kernel_list, &orbit_no,
                                     &orbit_length);

               /* edges in the graph */
               i__ = data->first_aff[i];
               if (orbit_no > 0){
                  gitter_no[0][0]++;
                  gitter_no[0][gitter_no[0][0]] = k;
               }

               for (m = 0; m < orbit_no; m++){
                  j__ = number_of_affine_class(data, orbit_rep[m], j, 0, oflag, &konj_wort);
                  zwischenerg[k].x[0][j__]++;
                  if (zwischenerg[k].l[0][j__] == NULL){
                     zwischenerg[k].l[0][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
                  }
		  if (oflag && zwischenerg[k].o[0][j__] == NULL){
                     zwischenerg[k].o[0][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
		  }
                  zwischenerg[k].l[0][j__][0]++;
                  zwischenerg[k].l[0][j__][ zwischenerg[k].l[0][j__][0] ] =
                     orbit_length[m];
		  if (oflag){
		     zwischenerg[k].o[0][j__][0]++;
		     zwischenerg[k].o[0][j__][ zwischenerg[k].o[0][j__][0] ] =
		        obergruppenzahl(invlattice, data->Z[j]->normal,
			                data->stab_coz[j][j__],
					data->stab_gen_no[j][j__], konj_wort);
                     free(konj_wort);
                  }
               }
               free(orbit_length);
               for (m = 0; m < orbit_no; m++){
                  free_mat(orbit_rep[m]);
               }
               free(orbit_rep);


               /* l > 0 (d.h. Untergruppen der nicht zerfallenden Gruppen) */
	       /* orbit on ksi + Ker(phi): l = 0 */
               kernel_elements_2_affine(kernel_elements, data->coho_size[j]);

               for (l = 1; l < reps_no; l++){
	          /* welcher Repraesentant ist Untergruppe des
		     Standardvertreter einer affinen Klasse (welcher?)? */
                  reps_flag = 0;
                  for (m = 0; m < length[l]; m++){
                     coz_i = mat_mul(phi, reps[l][m]);
                     i__ = number_of_affine_class(data, coz_i, i, 1, FALSE, NULL);
                     free_mat(coz_i);
                     if (m == 0)
                        i__first = i__;
                     if (abs(i__first) != abs(i__)){
                        /* paranoia test */
                        fprintf(stderr, "ERROR 2 in subgroupgraph!\n");
                        exit(9);
                     }
                     if (i__ > 0){
                        reps_flag = 1;
                        /* break; */
                     }
                  }

                  if (reps_flag != 1){
                     /* printf("HHHHHHHHHHHHHHHHHHHHH\n"); */
		     /* Die Gruppen U_{l,m}, die gegeben sind durch reps[l][m],
		        sind Untergruppen einer Raumgruppe R, die nicht der
			Standardvertreter S der affinen Klasse ist (fuer alle m)!
			Es sei n mit R^n = S. Dann ist U_l^n Untergruppe
			von S. Diese Gruppen werden bei dem Gitter nL betrachtet,
			wobei L das aktuell betrachtete Gitter sei. L und nL sind
			offensichtlich nicht in einer Bahn unter der Punktgruppe
			des affinen Nomalisators von S, d.h. unter dem Stabilisator
			des Coykels von S. */
		  }
		  else{
                     i__ = abs(i__); /* nur noetig, wenn break auskommentiert */

                     /* consider only standard representatives */
                     counter = 0;
                     S_ksi = (matrix_TYP **)calloc(WORDS_no[l], sizeof(matrix_TYP *));
                     for (m = 0; m < WORDS_no[l]; m++){
                        /* simplify the words */
                        normalize_word(WORDS[l][m]);
                        if (WORDS[l][m][0] == 1 && WORDS[l][m][1] < 0)
                           WORDS[l][m][1] *= (-1);
                        if (word_already_there(WORDS[l], m) == 0){
                           S_ksi[counter] = graph_mapped_word(WORDS[l][m], S1, S1_inv, data->X[j][1]);
                           if (mat_search(S_ksi[counter], S_ksi, counter, mat_comp) != -1){
                              free_mat(S_ksi[counter]);
                           }
                           else{
                              mat_quicksort(S_ksi, 0, counter, mat_comp);
                              counter++;
                           }
                        }
                     }

		     /* Da es uns nur um Anzahlen geht, koennen wir auch o.B.d.A. Bahnen auf
		        reps[l][0] + Kern (phi) berechnen */
                     orbit_rep = orbit_ksi_plus_ker(reps[l][0], kernel_elements, kernel_order,
                                                    data->X[j][1], S_ksi, counter,
                                                    data->coho_size[j], kernel_list, &orbit_no,
                                                    &orbit_length);

                     /* edges in tbe graph */
                     if (gitter_no[i__][gitter_no[i__][0]] != k || gitter_no[i__][0] == 0){
                        gitter_no[i__][0]++;
                        gitter_no[i__][gitter_no[i__][0]] = k;
                     }
                     for (m = 0; m < orbit_no; m++){
                        coz_j = mat_add(orbit_rep[m], reps[l][0], eins, eins);
                        j__ = number_of_affine_class(data, coz_j, j, 0, oflag, &konj_wort);
                        free_mat(coz_j);
                        zwischenerg[k].x[i__][j__]++;
                        if (zwischenerg[k].l[i__][j__] == NULL){
                           zwischenerg[k].l[i__][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
                        }
                        if (oflag && zwischenerg[k].o[i__][j__] == NULL){
                           zwischenerg[k].o[i__][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
                        }
                        zwischenerg[k].l[i__][j__][0]++;
                        zwischenerg[k].l[i__][j__][ zwischenerg[k].l[i__][j__][0] ] =
                           orbit_length[m];
			if (oflag){
                           zwischenerg[k].o[i__][j__][0]++;
                           zwischenerg[k].o[i__][j__][ zwischenerg[k].o[i__][j__][0] ] =
                              obergruppenzahl(invlattice, data->Z[j]->normal,
			                      data->stab_coz[j][j__],
					      data->stab_gen_no[j][j__], konj_wort);
			   free(konj_wort);
			}
                     }

                     /* clean */
                     free(orbit_length);
                     for (m = 0; m < orbit_no; m++){
                        free_mat(orbit_rep[m]);
                     }
                     free(orbit_rep);
                     for (m = 0; m < counter; m++)
                        free_mat(S_ksi[m]);
                     free(S_ksi);
                  }
               }
            }
            else{
               /* trivial part */
               switch (H1_mod_ker.flag){
                  case 0:
                  case 2:
                     /* orbit on 0 + Ker(phi) */
                     orbit_rep = orbit_ker(kernel_elements, kernel_order, data->X[j][1], S1,
                                           S1_word_no, data->coho_size[j], kernel_list, &orbit_no,
                                           &orbit_length);

                     /* edges in the graph */
                     i__ = data->first_aff[i];
                     if (orbit_no > 0){
                        gitter_no[0][0]++;
                        gitter_no[0][gitter_no[0][0]] = k;
                     }
                     for (m = 0; m < orbit_no; m++){
                        j__ = number_of_affine_class(data, orbit_rep[m], j, 0, oflag, &konj_wort);
                        zwischenerg[k].x[0][j__]++;
                        if (zwischenerg[k].l[0][j__] == NULL){
                           zwischenerg[k].l[0][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
                        }
                        if (oflag && zwischenerg[k].o[0][j__] == NULL){
                           zwischenerg[k].o[0][j__] = (int *)calloc(orbit_no + 1, sizeof(int));
                        }
                        zwischenerg[k].l[0][j__][0]++;
                        zwischenerg[k].l[0][j__][ zwischenerg[k].l[0][j__][0] ] =
                           orbit_length[m];
			if (oflag){
                           zwischenerg[k].o[0][j__][0]++;
                           zwischenerg[k].o[0][j__][ zwischenerg[k].o[0][j__][0] ] =
                              obergruppenzahl(invlattice, data->Z[j]->normal,
			                      data->stab_coz[j][j__],
					      data->stab_gen_no[j][j__], konj_wort);
			   free(konj_wort);
			}
                     }
                     free(orbit_length);
                     for (m = 0; m < orbit_no; m++){
                        free_mat(orbit_rep[m]);
                     }
                     free(orbit_rep);
                     break;
                  case 1:
                  case 3:
                     zwischenerg[k].x[0][0]++;
                     zwischenerg[k].l[0][0] = (int *)calloc(2, sizeof(int));
                     zwischenerg[k].l[0][0][0] = 1;
                     zwischenerg[k].l[0][0][1] = 1;
		     if (oflag){
                        zwischenerg[k].o[0][0] = (int *)calloc(2, sizeof(int));
                        zwischenerg[k].o[0][0][0] = 1;
                        zwischenerg[k].o[0][0][1] =
			   obergruppenzahl(invlattice, data->Z[j]->normal,
			                   data->stab_coz[j][0],
				   	   data->stab_gen_no[j][0], NULL);
		     }
                     gitter_no[0][0]++;
                     gitter_no[0][gitter_no[0][0]] = k;
                     break;
                  default:
                     fprintf(stderr, "ERROR 4 in subgroupgraph!\n");
                     exit(5);
               }
            }

            /* clean */
            if (H1_mod_ker.flag == 0 && H1_mod_ker.erz_no > 0){
               for (l = 0; l < S1_word_no; l++){
                  free_mat(new_rep[l]);
               }
               free(new_rep);
               for (l = 1; l < reps_no; l++){
                  for (m = 0; m < WORDS_no[l]; m++){
                     free(WORDS[l][m]);
                  }
                  for (m = 0; m < length[l]; m++){
                     free_mat(reps[l][m]);
                  }
                  free(reps[l]);
                  free(WORDS[l]);
               }
               free(WORDS);
               free(WORDS_no);
               free(reps);
               free(length);
            }
            free_H1_mod_ker_TYP(H1_mod_ker);
            for (l = 0; l < data->coho_size[j]; l++){
               if (kernel_elements[l] != NULL)
                  free_mat(kernel_elements[l]);
            }
            free(kernel_elements);
            free(kernel_list);
            for (l = 0; l < kernel_mat->cols; l++){
               free_mat(kernel_gen[l]);
            }
            free(kernel_gen);
            free_mat(kernel_mat);
            for (l = 0; l < image_gen_no; l++){
               free_mat(image[l]);
            }
            if (phi->cols > 0){
               for (l = 0; l < S1_word_no; l++){
                  free_mat(S1[l]);
                  if (S1_inv[l] != NULL)
                     free_mat(S1_inv[l]);
               }
               free(S1);
               free(S1_inv);
            }
            free(image);
            free_mat(diag);
            free_mat(invlattice);
            free_mat(phi);
         }

         /* output */
         for (k = 0; k < data->aff_no[i] && entry.anz > 1; k++){

	   /* mehrere Teilgitter */
            if (gitter_no[k][0] > 1){
               /* orbit on lattices for each affine class with subgroups from
                  different lattices */
               lattice_orbits = (int *)calloc(entry.anz, sizeof(int));
               lattice_orbit_length = (int *)calloc(entry.anz + 1, sizeof(int));
               smallest = (int *)calloc(entry.anz + 1, sizeof(int));
               /* Stabilisatoren nur bei Bedarf berechnen (get_Q_data und
                  free_Q_data aendern !!!!!)
               if (data->stab_coz[i][k] == NULL){
                  data->stab_coz[i][k] = stab_coz(data->WORDS[i][k], data->NUMBER_OF_WORDS[i][k],
                                                  data->Z[i]->normal, data->norm_inv[i],
                                                  data->Z[i]->normal_no, j,
                                                  &data->stab_gen_no[i][k]);
               }
               */
               lattice_orbits_no = orbit_on_lattices(entry.lsf, entry.anz, data->stab_coz[i][k],
                                                     data->stab_gen_no[i][k], lattice_orbits,
                                                     lattice_orbit_length, smallest, NULL);

               for (l = 1; l <= lattice_orbits_no; l++){
                  if (gibt_gitter_beitrag(gitter_no[k], gitter_no[k][0], smallest[l])){
                     printf("%i: ", k + 1 + data->first_aff[i]);
                     for (m = 0; m < data->aff_no[j]; m++){
                        if (zwischenerg[smallest[l]].x[k][m] > 0){
                           for (q = 0; q < zwischenerg[smallest[l]].x[k][m]; q++){
                              printf("%i (%i, ", m + 1 + data->first_aff[j],
                                                lattice_orbit_length[l] *
                                                zwischenerg[smallest[l]].l[k][m][q + 1] *
                                                kernel_factor[smallest[l]]);
			      if (oflag)
			         printf("%i, ", zwischenerg[smallest[l]].o[k][m][q + 1]);
                              printf("%i^%i)  ", primes[smallest[l]], exponent[smallest[l]]);
                              erg->array.SZ[k + data->first_aff[i]][m + data->first_aff[j]]++;
                           }
                        }
                     }
                     printf("\n");
                  }
               }

               free(lattice_orbits);
               free(lattice_orbit_length);
               free(smallest);
            }
            else{
               /* for this affine class, there are only subgroups from one lattice */
               if (gitter_no[k][0] == 1){
                  printf("%i: ", k + 1 + data->first_aff[i]);
                  for (l = 0; l < data->aff_no[j]; l++){
                     if (zwischenerg[gitter_no[k][1]].x[k][l] > 0){
                        for (q = 0; q < zwischenerg[gitter_no[k][1]].x[k][l]; q++){
                           printf("%i (%i, ", l + 1 + data->first_aff[j],
                                               zwischenerg[gitter_no[k][1]].l[k][l][q + 1] *
                                               kernel_factor[gitter_no[k][1]]);
			   if (oflag)
                              printf("%i, ", zwischenerg[gitter_no[k][1]].o[k][l][q + 1]);
                           printf("%i^%i)  ", primes[gitter_no[k][1]], exponent[gitter_no[k][1]]);
                           erg->array.SZ[k + data->first_aff[i]][l + data->first_aff[j]]++;
                        }
                     }
                  }
                  printf("\n");
               }
            }
         }
         if (entry.anz == 1){
            /* there is only one lattice */
            for (k = 0; k < data->aff_no[i]; k++){
               if (gitter_no[k][0] != 0){
                  printf("%i: ", k + 1 + data->first_aff[i]);
                  for (l = 0; l < data->aff_no[j]; l++){
                     if (zwischenerg[0].x[k][l] > 0){
                        for (q = 0; q < zwischenerg[0].x[k][l]; q++){
			   printf("%i (%i, ", l + 1 + data->first_aff[j],
			                      zwischenerg[0].l[k][l][q + 1] *
                                              kernel_factor[0]);
                           if (oflag)
			       printf("%i, ", zwischenerg[0].o[k][l][q + 1]);
                           printf("%i^%i)  ", primes[0], exponent[0]);
                           erg->array.SZ[k + data->first_aff[i]][l + data->first_aff[j]]++;
                        }
                     }
                  }
                  printf("\n");
               }
            }
         }
         free_SOL_TYP(zwischenerg, entry.anz, data->aff_no[i], data->aff_no[j]);
         for (k = 0 ; k < data->aff_no[i]; k++){
            free(gitter_no[k]);
         }
         free(gitter_no);
         free(kernel_factor);
         free(exponent);
         free(primes);
      }
      for (j = 0; j < translanz; j++){
         free_mat(translationen[j]);
      }
      free(translationen);
   }

   /* clean up */
   for (i = 0; i < data->pres->rows; i++)
      wordfree(relator + i);
   free(relator);
   free_mat(id);
   for (i = 0; i < data->Z_no; i++){
      for (j = 0; j < data->Z[i]->normal_no; j++){
         if (dataNinv[i][j] != NULL)
            free_mat(dataNinv[i][j]);
      }
      free(dataNinv[i]);
   }
   free(dataNinv);

   return(erg);
}




