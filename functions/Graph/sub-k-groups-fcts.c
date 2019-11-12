/* last change: 07.02.2001 by Oliver Heidbuechel */


#include "ZZ.h"
#include "typedef.h"
#include "utils.h"
#include "presentation.h"
#include "matrix.h"
#include "bravais.h"
#include "base.h"
#include "datei.h"
#include "graph.h"
#include "gmp.h"
#include "zass.h"
#include "tools.h"
#include "longtools.h"
#include "voronoi.h"
#include "symm.h"
#include "autgrp.h"
#include "reduction.h"

#define MYSIZE 1024


/* --------------------------------------------------------------------- */
/* Konjugiere bravais_TYP                                                */
/* (genauer berechne T B T^{-1})                                         */
/* Quellcode konj_bravais.c und normalizer.c entnommen                   */
/* Es wird nur ->gen, ->gen_no, ->normal, ->normal_no, ->order und       */
/* ->divisors, ->form, ->form_no angelegt.                               */
/* --------------------------------------------------------------------- */
/* B: zu konjugierende Gruppe                                            */
/* T: Matrix mit der konjugiert wird                                     */
/* Ti: Inverse zu T                                                      */
/* --------------------------------------------------------------------- */
static bravais_TYP *my_konj_bravais(bravais_TYP *B,
                                    matrix_TYP *T,
				    matrix_TYP *Ti)
{
   int i, Vanz;
   bravais_TYP *G, *BG, *BGtr;
   matrix_TYP *waste, *tmp, *tmp2, *A;
   voronoi_TYP **V;

   G = init_bravais(B->dim);

   /* Order */
   if (B->order != 0){
     G->order = B->order;
     memcpy(G->divisors,B->divisors,100*sizeof(int));
   }

   /* Generators */
   if(B->gen_no != 0)
   {
     G->gen_no = B->gen_no;
     G->gen = (matrix_TYP **)xmalloc(B->gen_no *sizeof(matrix_TYP *));
     for(i=0;i<B->gen_no;i++)
     {
        waste = mat_mul(T, B->gen[i]);
        G->gen[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
   }

   /* Normalizer, weiter wird der Formenraum von G angelegt */
   BG = bravais_group(G, FALSE);

   /* let's see whether we already got the formspace */
   if (BG->form == NULL){
      BG->form = formspace(BG->gen, BG->gen_no, 1, &BG->form_no);
   }
   BGtr = tr_bravais(BG, 1, FALSE);

   /* firstly calculate an positive definite BG-invariant form */
   tmp2 = init_mat(BG->dim,BG->dim,"1");
   tmp = rform(BG->gen,BG->gen_no,tmp2,101);
   free_mat(tmp2);

   /* now calculate the trace bifo */
   tmp2 = trace_bifo(BG->form,BGtr->form,BG->form_no);
   A = first_perfect(tmp, BG, BGtr->form, tmp2, &Vanz);
   free_mat(tmp2);
   free_mat(tmp);

   V = normalizer(A, BG, BGtr, 1949, &Vanz);

   /* now we got BG and it's normalizer, so see what we can do with G */
   G->normal = normalizer_in_N(G, BG, &G->normal_no, FALSE);

   /* clean */
   for(i = 0; i < Vanz; i++){
      clear_voronoi(V[i]);
      free(V[i]);
   }
   free(V);
   free_mat(A);
   free_bravais(BG);
   free_bravais(BGtr);

   return(G);
}





/* --------------------------------------------------------------------- */
/* Calculate H^1(G, Q^n/Z^n)                                             */
/* G < GL_n(Z) finite, relator: presentation in relator format           */
/* relatoranz: number of rows in the presentation                        */
/* see .../Zassen/zass.c                                                 */
/* --------------------------------------------------------------------- */
matrix_TYP **calculate_H1(bravais_TYP *G,
                          word *relator,
                          int relatoranz)
{
   matrix_TYP **X,
              **matinv;

   long dim;

   int i;


   /* prepare */
   matinv = (matrix_TYP **)calloc(G->gen_no, sizeof(matrix_TYP *));

   /* calculate H^1 */
   X = cohomology(&dim, G->gen, matinv, relator, G->gen_no, relatoranz);

   /* clean */
   for (i = 0; i < G->gen_no; i++){
      if (matinv[i] != NULL)
         free_mat(matinv[i]);
   }
   free(matinv);

   return(X);
}



/* --------------------------------------------------------------------- */
/* free the structure coz_TYP                                            */
/* --------------------------------------------------------------------- */
void free_coz_TYP(coz_TYP coz_info)
{
   int i;

   mpz_clear(&coz_info.name);
   mpz_clear(&coz_info.number);
   free_mat(coz_info.std_coz);
   free_mat(coz_info.diff);
   for (i = 0; coz_info.WORDS != NULL && i < coz_info.WORDS_no; i++)
      free(coz_info.WORDS[i]);
   if (coz_info.WORDS != NULL)
      free(coz_info.WORDS);
   if (coz_info.Stab){
      for (i = 0; i < coz_info.Stab_no; i++)
         free_mat(coz_info.Stab[i]);
      free(coz_info.Stab);
   }
   if (coz_info.darst != NULL)
      free_mat(coz_info.darst);
   for (i = 0; coz_info.aff_transl != NULL && i < coz_info.aff_transl_no; i++)
      free_mat(coz_info.aff_transl[i]);
   if (coz_info.aff_transl)
      free(coz_info.aff_transl);
   if (coz_info.N_darst){
      for (i = 0; i < coz_info.N_darst_no; i++)
         free_mat(coz_info.N_darst[i]);
      free(coz_info.N_darst);
   }
   for (i = 0; i < coz_info.N_darst_no; i++){
      if (coz_info.N_inv[i] != NULL) free_mat(coz_info.N_inv[i]);
   }
   free(coz_info.N_inv);
   free_mat(coz_info.coz);
   for (i = 0; i < 3; i++)
      free_mat(coz_info.X[i]);
   free(coz_info.X);
   free_mat(coz_info.GLS);
}



/* --------------------------------------------------------------------- */
/* Calculate information about a cocycle, i. e. fill the structure       */
/* coz_TYP for a coz given in coz, G is the corresponding point-group    */
/* X contains the return value of cohomology for G                       */
/* G->normal has to generate together with G->gen N_Gl_n(Z) (G).         */
/* G->zentr isn't considered!!!                                          */
/* --------------------------------------------------------------------- */
/* coz_info.Stab = Stabilisator des Cozykels / G !!!                     */
/* --------------------------------------------------------------------- */
coz_TYP identify_coz(bravais_TYP *G,
                     matrix_TYP *coz,
                     matrix_TYP **X)
{
   matrix_TYP **TR, **N_inv;

   coz_TYP coz_info;

   rational eins, minuseins;

   MP_INT null;

   int i, zen_no;


   zen_no = G->zentr_no;
   G->zentr_no = 0;
   coz_info.GLS = mat_inv(X[2]);

   if (X[0]->cols >= 1){
      /* prepare */
      eins.z = eins.n = minuseins.n = 1;
      minuseins.z = -1;
      mpz_init_set_si(&coz_info.name, 0);
      mpz_init_set_si(&coz_info.number, 0);
      mpz_init_set_si(&null, 0);
      coz_info.WORDS = (int **)calloc(MIN_SPEICHER, sizeof(int *));
      coz_info.WORDS_no = 0;
      N_inv = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));

      /* get information about the cocycle */
      TR = identify(X[0], X[1], X[2], G, &coz, &coz_info.name, 1, 1,
                    &coz_info.WORDS, &coz_info.WORDS_no);

      coz_info.darst = standard_rep(coz, coz_info.GLS, X[1]);
      valuation(coz_info.darst, X[1], &coz_info.number);
      coz_info.std_coz = convert_to_cocycle(coz_info.darst, X[0], X[1]);
      coz_info.diff = mat_add(coz, coz_info.std_coz, eins, minuseins);

      coz_info.Stab = stab_coz(coz_info.WORDS, coz_info.WORDS_no, G->normal,
                               N_inv, G->normal_no, mpz_cmp(&coz_info.name, &null),
                               &coz_info.Stab_no);

      /* represenation of G->normal on H^1(G,Q^n/Z^n) */
      coz_info.N_darst_no = G->normal_no;
      coz_info.N_darst = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));

      /* eigentlich UEBERFLUESSIG:wird schon in identify berechnet */ 
      for (i = 0; i < G->normal_no; i++){
         coz_info.N_darst[i] = normalop(X[0], X[1], X[2], G, G->normal[i], 1);
      }

      /* clean */
      free_mat(TR[0]);
      free(TR);
      mpz_clear(&null);
      for (i = 0; i < G->normal_no; i++){
         if (N_inv[i] != NULL) free_mat(N_inv[i]);
      }
      free(N_inv);
   }
   else{
      /* trivial H^1 */
      coz_info.std_coz = init_mat(coz->rows, 1, "");
      mpz_init_set_si(&coz_info.name, 0);
      mpz_init_set_si(&coz_info.number, 0);
      coz_info.diff = copy_mat(coz);
      coz_info.WORDS = NULL;
      coz_info.Stab = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));
      for (i = 0; i < G->normal_no; i++){
         coz_info.Stab[i] = copy_mat(G->normal[i]);
      }
      coz_info.Stab_no = G->normal_no;
      coz_info.darst = init_mat(0, 0, "");
      coz_info.N_darst = NULL;
      coz_info.N_darst_no = 0;
   }
   G->zentr_no = zen_no;
   coz_info.aff_transl = NULL;
   coz_info.aff_transl_no = 0;
   coz_info.coz = copy_mat(coz);
   coz_info.X = (matrix_TYP **)calloc(3, sizeof(matrix_TYP *));
   coz_info.N_inv = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));
   for (i = 0; i < 3; i++)
      coz_info.X[i] = copy_mat(X[i]);

   return(coz_info);
}



/* --------------------------------------------------------------------- */
/* informations about G-invariant sublattices                            */
/* --------------------------------------------------------------------- */
typedef struct
{
   matrix_TYP **lattices;	/* sublattices (matrices in long_col_hnf - Form) */
   boolean *trivialflag;	/* trivialflag == TRUE iff lattice[i] == p_i * Z^n */
   int no;                      /* number of sublattices */
   int *orbitlist;		/* lattices[i] belongs to orbit orbitlist[i] in 1,2, ... */
   int *orbitlength;		/* orbitlength[i]: length of the orbits i = 1, 2, ... */
   int *smallest;		/* smallest[i]: number of the lattice with the lowest
                                   number in orbit[i], i = 1, 2, ... */
   int orbit_no;		/* number of orbits */
   matrix_TYP **conj;		/* conj[i] * lattices[smallest[orbitlist[i]]] = lattices[i] */
} sublattice_TYP;



/* --------------------------------------------------------------------- */
/* free the structure sublattice_TYP                                     */
/* --------------------------------------------------------------------- */
static void free_sublattice_TYP(sublattice_TYP SL)
{
   int i;


   for (i = 0; i < SL.no; i++){
      free_mat(SL.lattices[i]);
      free_mat(SL.conj[i]);
   }
   free(SL.lattices);
   free(SL.conj);
   free(SL.orbitlist);
   free(SL.orbitlength);
   free(SL.smallest);
   free(SL.trivialflag);
}



/* --------------------------------------------------------------------- */
/* returns NULL iff cocycle "coz" isn't in the image or is 0             */
/* otherwise returns sol with phi(sol) = coz_darst                       */
/* (C_a x C_b x ... - representation)                                    */
/* --------------------------------------------------------------------- */
static matrix_TYP *cocycle_in_image(matrix_TYP *coz_darst,
                                    matrix_TYP **image_gen,
                                    int image_gen_no,
                                    matrix_TYP *D,
                                    int flag,
                                    int erz_no,
                                    boolean *is_in_image)
{
   matrix_TYP *sol = NULL,
              **images, **koords,
              *tmp;

   int i, j, k, pos, first, last, diff,
       anz, size = MYSIZE,
       *number;

   MP_INT val;



   /* two trivial cases */
   if (flag == 2 || flag == 3 || equal_zero(coz_darst)){
      is_in_image[0] = TRUE;
      return(sol);
   }

   /* further cases */
   if (flag == 1){
      is_in_image[0]= FALSE;
   }
   else if (flag == 0){
      if (erz_no == 0){
         is_in_image[0] = FALSE;
      }
      else{
         /* prepare */
         anz = 1;
         for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);
         for (last = first; last < D->cols && D->array.SZ[last][last] != 0; last++);
         diff = last - first;
         images = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
         images[0] = init_mat(diff, 1, "");
         number = (int *)calloc(size, sizeof(int));
         koords = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
         koords[0] = init_mat(image_gen_no, 1, "");
         mpz_init(&val);
         is_in_image[0] = FALSE;

         /* orbit algorithm (surely, there are better methods)*/
         for (i = 0; !is_in_image[0] && i < anz; i++){
            for (j = 0; !is_in_image[0] && j < image_gen_no; j++){
               if (!equal_zero(image_gen[j])){
                  tmp = add_mod_D(images[i], image_gen[j], D, first, diff);
                  valuation(tmp, D, &val);
                  pos = mpz_get_ui(&val);
                  for (k = 0; k < anz; k++){
                     if (number[k] == pos)
                       break;
                  }
                  if (k == anz){
                     /* new element in the orbit */
                     if (anz >= size){
                        size += MYSIZE;
                        images = (matrix_TYP **)realloc(images, size * sizeof(matrix_TYP *));
                        number = (int *)realloc(number, size * sizeof(int));
                        koords = (matrix_TYP **)realloc(koords, size * sizeof(matrix_TYP *));
                     }
                     images[anz] = tmp;
                     koords[anz] = copy_mat(koords[i]);
                     koords[anz]->array.SZ[j][0]++;
                     number[anz] = pos;
                     if (cmp_mat(images[anz], coz_darst) == 0){
                        sol = copy_mat(koords[anz]);
                        is_in_image[0] = TRUE;
                     }
                     anz++;
                  }
                  else
                     free_mat(tmp);
               }
            }
         }

         /* clean */
         for (i = 0; i < anz; i++){
            free_mat(koords[i]);
            free_mat(images[i]);
         }
         free(koords);
         free(images);
         free(number);
         mpz_clear(&val);
      }
   }
   else{
      fprintf(stderr, "This flag isn't allowed. ERROR in cocycle_in_image!\n");
      exit(4);
   }

   return(sol);
}



/* --------------------------------------------------------------------- */
/* H1: return value of cohomology                                        */
/* koord: koordinates of a cocyle in H^1(...) = C_a x C_b x ...          */
/* returns cocycle                                                       */
/* --------------------------------------------------------------------- */
matrix_TYP *darst_to_C1(matrix_TYP **H1,
                        matrix_TYP *koord)
{
   matrix_TYP *C1, *tmp;

   int first, last, diff, i, j;

   rational eins, zahl;


   /* prepare */
   for (first = 0; first < H1[1]->cols && H1[1]->array.SZ[first][first] == 1; first++);
   for (last = first; last < H1[1]->cols && H1[1]->array.SZ[last][last] != 0; last++);
   diff = last - first;
   eins.z = eins.n = 1;

   /* paranoia */
   if (diff != koord->rows){
      fprintf(stderr, "ERROR in darst_to_C1!!!\n");
      exit(9);
   }

   C1 = init_mat(H1[0]->rows, 1, "");
   for (i = 0; i < diff; i++){
      if (koord->array.SZ[i][0] != 0){
         tmp = init_mat(H1[0]->rows, 1, "");
         for (j = 0; j < H1[0]->rows; j++){
            tmp->array.SZ[j][0] = H1[0]->array.SZ[j][i];
         }
         zahl.z = koord->array.SZ[i][0];
         zahl.n = H1[1]->array.SZ[i + first][i + first];
         Normal(&zahl);
         C1 = mat_addeq(C1, tmp, eins, zahl);
         free_mat(tmp);
      }
   }

   return(C1);
}



/* --------------------------------------------------------------------- */
/* adds an element of B^1(G, Q^n/L) to preimage such that afterwarts     */
/* phi (preimage) = coz (in C^1)                                         */
/* case: normalizer element is ID                                        */
/* [part of this code is part of the code in coboundary]                 */
/* --------------------------------------------------------------------- */
static void korrektes_urbild_ID(matrix_TYP *preimage,
                                matrix_TYP *coz,
                                bravais_TYP *G)
{
   matrix_TYP *diff, *A, *B, *tmp, **erg;

   rational eins, minuseins;

   int i, n, k, l;


   eins.z = eins.n = minuseins.n = 1; minuseins.z = -1;
   diff = mat_add(coz, preimage, eins, minuseins);

   n = G->dim;

   B = init_mat(n * G->gen_no, n, "k");

   /* belegen der matrix B = (g_1 - I, ...)^tr */
   for (i = 0; i < G->gen_no; i++){
      for (k = 0; k < n; k++){
         for (l = 0; l < n; l++){
    	    if (k == l){
               B->array.SZ[k+i*n][l] = G->gen[i]->array.SZ[k][l] - 1;
            }
            else{
               B->array.SZ[k+i*n][l] = G->gen[i]->array.SZ[k][l];
            }
         }
      }
   }

   A = copy_mat(B);
   real_mat(A, A->rows, A->cols + 1);

   /* rechte Seite */
   for (i = 0; i < G->gen_no; i++)
      for (k = 0; k < n; k++)
         A->array.SZ[i * n + k][n] = diff->array.SZ[i * n + k][0];

   k = long_row_gauss(A);
   real_mat(A, k, A->cols);

   /* kick out the rows which have 0's in the first n entries */
   for (i = A->rows - 1; i > 0; i--){
      for (k = 0; k < n && A->array.SZ[i][k] == 0; k++);
      if (k == n){
         if (A->array.SZ[i][n] % diff->kgv){
            fprintf(stderr,"error in korrektes_urbild_ID\n");
            exit(3);
         }
         else{
            real_mat(A,A->rows-1,A->cols);
         }
      }
   }

   real_mat(diff, A->rows, 1);

   for (i = 0; i < A->rows; i++)
      diff->array.SZ[i][0] = A->array.SZ[i][n];

   real_mat(A,A->rows, n);

   Check_mat(A);
   Check_mat(diff);

   /* put_mat(A,0,0,0);
   put_mat(diff,0,0,0); */

   erg = long_solve_mat(A, diff);

   if (erg == NULL || erg[0] == NULL){
      fprintf(stderr,"error in korrektes_urbild_ID\n");
      exit(3);
   }

   tmp = mat_mul(B, erg[0]);

   mat_addeq(preimage, tmp, eins, eins);

   /* clean */
   if (erg[0] != NULL) free_mat(erg[0]);
   if (erg[1] != NULL) free_mat(erg[1]);
   free(erg);
   free_mat(A);
   free_mat(B);
   free_mat(diff);
   free_mat(tmp);
}


/* --------------------------------------------------------------------- */
/* adds an element of B^1(G, Q^n/L) to preimage such that afterwarts     */
/* phi (preimage) = coz (in C^1)                                         */
/* case: normalizer element is possibly NOT ID                           */
/* [part of this code is part of the code in coboundary]                 */
/* --------------------------------------------------------------------- */
static void korrektes_urbild(matrix_TYP **preimage,
                             coz_TYP *coz_info,
                             bravais_TYP *G,
                             matrix_TYP *L,
                             boolean flag)
{
   matrix_TYP *rep, *paranoia, **N_darst_inv, **N, **NNN, **orbit, *n, *tmp, **conj;

   bravais_TYP *R, *RR;

   int **words, word_no, i, j, k,
       anz, size = MYSIZE, first, last, diff, pos,
       *number;

   boolean found;

   MP_INT val;



   rep = standard_rep(preimage[0], coz_info->GLS, coz_info->X[1]);

   if (cmp_mat(rep, coz_info->darst) == 0){
      korrektes_urbild_ID(preimage[0], coz_info->coz, G);
      free_mat(rep);
   }
   else{

      /* calculate the stabilizer of the lattice */
      if (flag){
         N = coz_info->N_darst;
         NNN = G->normal;
         word_no = G->normal_no;
      }
      else{
         words = stab_lattice(L, G->normal, G->normal_no, &word_no, NULL, 0, NULL);

         N_darst_inv = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));
         N = (matrix_TYP **)calloc(word_no, sizeof(matrix_TYP *));
         NNN = (matrix_TYP **)calloc(word_no, sizeof(matrix_TYP *));
         for (i = 0; i < word_no; i++){
            NNN[i] = mapped_word(words[i], G->normal, coz_info->N_inv);
            paranoia = mat_mul(NNN[i], L);
            long_col_hnf(paranoia);
            if (cmp_mat(paranoia, L) != 0){
               fprintf(stderr, "ERROR in korrektes_urbild!\n");
               exit(9);
            }
            free_mat(paranoia);
            N[i] = graph_mapped_word(words[i], coz_info->N_darst, N_darst_inv, coz_info->X[1]);
            free(words[i]);
         }
         free(words);
         for (i = 0; i < G->normal_no; i++){
            if (N_darst_inv[i] != NULL) free_mat(N_darst_inv[i]);
         }
         free(N_darst_inv);
      }

      /* find n in N with n rep = coz_info->darst */
      for (first = 0; first < coz_info->X[1]->cols && coz_info->X[1]->array.SZ[first][first] == 1; first++);
      for (last = first; last < coz_info->X[1]->cols && coz_info->X[1]->array.SZ[last][last] != 0; last++);
      diff = last - first;
      found = FALSE;
      anz = 1;
      orbit = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
      number = (int *)calloc(size, sizeof(int));
      conj = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
      conj[0] = init_mat(G->dim, G->dim, "1");
      orbit[0] = copy_mat(rep);
      tmp = NULL;
      mpz_init(&val);
      valuation(orbit[0], coz_info->X[1], &val);
      number[0] = mpz_get_ui(&val);
      n = NULL;

      /* out of orbit_rep */
      for (i = 0; !found && i < anz; i++){
         for (j = 0; !found && j < word_no; j++){
            tmp = mat_mul(N[j], orbit[i]);
            for (k = 0; k < diff; k++){
               tmp->array.SZ[k][0] %= coz_info->X[1]->array.SZ[k + first][k + first];
               if (tmp->array.SZ[k][0] < 0)
                  tmp->array.SZ[k][0] += coz_info->X[1]->array.SZ[k + first][k + first];
            }
            valuation(tmp, coz_info->X[1], &val);
            pos = mpz_get_ui(&val);
            for (k = 0; k < anz; k++){
               if (number[k] == pos)
               break;
            }
            if (k == anz){
               /* new element in the orbit */
               if (anz >= size){
                  size += MYSIZE;
                  orbit = (matrix_TYP **)realloc(orbit, size * sizeof(matrix_TYP *));
                  number = (int *)realloc(number, size * sizeof(int));
                  conj = (matrix_TYP **)realloc(conj, size * sizeof(matrix_TYP *));
               }
               orbit[anz] = tmp;
               conj[anz] = mat_mul(NNN[j], conj[i]);
               number[anz] = pos;
               if (cmp_mat(orbit[anz], coz_info->darst) == 0){
                  n = copy_mat(conj[anz]);
                  found = TRUE;
               }
               anz++;
            }
            else
               free_mat(tmp);
         }
      }

      /* paranoia test */
      if (n == NULL){
         fprintf(stderr, "ERROR 1 in korrektes_urbild\n");
         exit(3);
      }

      /* calculate the correct preimage */
      real_mat(n, n->rows + 1, n->cols);
      real_mat(n, n->rows, n->cols + 1);
      n->array.SZ[G->dim][G->dim] = 1;
      my_translation(n, preimage[0], coz_info->coz, G);
      R = extract_r(G, preimage[0]);
      free_mat(preimage[0]);
      RR = konj_bravais(R, n);
      preimage[0] = sg(RR, G);

      /* clean */
      free_mat(n);
      free_bravais(R);
      free_bravais(RR);
      for (i = 0; i < anz; i++){
         free_mat(conj[i]);
         free_mat(orbit[i]);
      }
      free(orbit);
      free(conj);
      free(number);
      mpz_clear(&val);
      if (flag){
         N = NULL;
         NNN = NULL;
      }
      else{
         for (i = 0; i < word_no; i++){
            free_mat(N[i]);
            free_mat(NNN[i]);
         }
         free(NNN);
         free(N);
      }
      free_mat(rep);
   }
}



/* --------------------------------------------------------------------- */
/* With one preimage and the kernels we can calculate all preimages in   */
/* C^1(G, Q^n/L)                                                         */
/* --------------------------------------------------------------------- */
static matrix_TYP **calculate_all_preimages(matrix_TYP *preimage,
                                            matrix_TYP **kernel_gen,
                                            matrix_TYP **k_g_darst,
                                            int kernel_gen_no,
                                            matrix_TYP *D,
                                            matrix_TYP **b1_kernel,
                                            int b1_kernel_no,
                                            int *anz)
{
   matrix_TYP **preimages, *tmp, **elements;

   int size = MYSIZE,
       i, j, k, pos, alt,
       *number,
       first, last, diff;

   MP_INT val;

   rational eins;


   eins.z = eins.n = 1;

   /* all preimages up to addition with element of Ker phi | B^1(B, Q^n/L) */
   if (kernel_gen_no > 0){
      /* prepare */
      anz[0] = 1;
      preimages = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
      preimages[0] = copy_mat(preimage);
      number = (int *)calloc(size, sizeof(int));
      elements = (matrix_TYP **)calloc(size, sizeof(matrix_TYP *));
      elements[0] = init_mat(k_g_darst[0]->rows, 1, "");
      mpz_init(&val);
      for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);
      for (last = first; last < D->cols && D->array.SZ[last][last] != 0; last++);
      diff = last - first;

      /* orbit algorithm */
      for (i = 0; i < anz[0]; i++){
         for (j = 0; j < kernel_gen_no; j++){
            tmp = add_mod_D(elements[i], k_g_darst[j], D, first, diff);
            valuation(tmp, D, &val);
            pos = mpz_get_ui(&val);
            for (k = 0; k < anz[0]; k++){
               if (number[k] == pos)
                  break;
            }
            if (k == anz[0]){
               /* new element in the orbit */
               anz[0]++;
               if (anz[0] >= size){
                  size += MYSIZE;
                  preimages = (matrix_TYP **)realloc(preimages, size * sizeof(matrix_TYP *));
                  number = (int *)realloc(number, size * sizeof(int));
                  elements = (matrix_TYP **)realloc(elements, size * sizeof(matrix_TYP *));
               }
               elements[anz[0] - 1] = tmp;
               preimages[anz[0] - 1] = mat_add(preimages[i], kernel_gen[j], eins, eins);
               number[anz[0] - 1] = pos;
            }
            else
               free_mat(tmp);
         }
      }

      /* clean */
      mpz_clear(&val);
      free(number);
      for (i = 0; i < anz[0]; i++)
         free_mat(elements[i]);
      free(elements);

      alt = anz[0];
      anz[0] = anz[0] * b1_kernel_no;
      if (size < anz[0])
         preimages = (matrix_TYP **)realloc(preimages, anz[0] * sizeof(matrix_TYP *));
   }
   else{
      /* trivial case */
      alt = 1;
      anz[0] = b1_kernel_no;
      preimages = (matrix_TYP **)calloc(anz[0], sizeof(matrix_TYP *));
      preimages[0] = copy_mat(preimage);
   }

   /* add elements of Ker phi | B^1(B, Q^n/L) */
   for (i = 1; i < b1_kernel_no; i++){
      for (j = 0; j < alt; j++){
         preimages[i * alt + j] = mat_add(preimages[j], b1_kernel[i], eins, eins);
      }
   }

   return(preimages);
}



/* --------------------------------------------------------------------- */
/* calculate the maximal klassengleich subgroups for R(G,coz) with       */
/* translation lattice L                                                 */
/* --------------------------------------------------------------------- */
/* G: Pointgroup, G->normal and G->gen have to generate N_{Gl_n(Z)}(G)   */
/* GL: G mit Gitter L konjugiert                                         */
/* aflag: calculate subgroups up to conjugation of the affine normalizer */
/*        of R                                                           */
/* anz: Speichere die Anzahl der Raumgruppen hier                        */
/* length: wenn aflag, dann speichere die Laengen der Bahnen hier        */
/* --------------------------------------------------------------------- */
static bravais_TYP **subgroups_for_L(coz_TYP *coz_info,
                                     matrix_TYP *L,
                                     bravais_TYP *G,
				     bravais_TYP *GL,
                                     matrix_TYP **H_G_Z,
                                     matrix_TYP **H_GL_Z,
                                     int *anz,
				     int **length,
				     boolean aflag)
{
   boolean is_in_image;

   int i, b1_kernel_no;

   matrix_TYP *diag, *GLS, *tmp, *tmp2, **reps,
              *phi, *kernel_mat, **kernel_gen, **image,
              *koord, *preimage, **all_preimages,
              **k, *nullcoz, **b1_kernel;

   H1_mod_ker_TYP H1_mod_ker;

   int image_gen_no;

   bravais_TYP **SG;


   /* prepare */
   anz[0] = 0;
   GLS = mat_inv(H_G_Z[2]);
   diag = matrix_on_diagonal(L, G->gen_no);

   /* calculate Phi */
   calculate_phi(diag, H_GL_Z[0], H_G_Z, H_GL_Z, GLS, &phi,
                 &kernel_mat, &image, &image_gen_no, &H1_mod_ker);


   /* is cocycle in the image (phi koord = coz->darst) */
   koord = cocycle_in_image(coz_info->darst, image, image_gen_no,
                            H_G_Z[1], H1_mod_ker.flag, H1_mod_ker.erz_no, &is_in_image);


   if (is_in_image){

      /* calculate information about Ker phi | B^1 */
      if (coz_info->aff_transl == NULL){ /* we calculate it only once */
         coz_info->aff_transl = transl_aff_normal(G->gen, G->gen_no, &coz_info->aff_transl_no);
      }
      b1_kernel = kernel_factor_fct(coz_info->aff_transl, coz_info->aff_transl_no, G->gen_no, L,
                                       &b1_kernel_no);

      if (!aflag){
         /* Berechne alle Untergruppen */
	 /* ========================== */

         /* calculate preimage of cocycle (in C^1)*/
         if (koord == NULL){
            preimage = init_mat(coz_info->coz->rows, 1, "");
         }
         else{
            tmp = darst_to_C1(H_GL_Z, koord);
            preimage = mat_mul(diag, tmp);
            free_mat(tmp);
         }
         korrektes_urbild(&preimage, coz_info, G, L, FALSE);

         /* calculate generators for Ker(phi on H^1(G,Q^n/L)) such that
            phi(gen) = 0 in C^1 */
         nullcoz = init_mat(coz_info->coz->rows, 1, "");
         kernel_gen = (matrix_TYP **)calloc(kernel_mat->cols, sizeof(matrix_TYP *));
         k = col_to_list(kernel_mat);
         for (i = 0; i < kernel_mat->cols; i++){
            tmp = darst_to_C1(H_GL_Z, k[i]);
            kernel_gen[i] = mat_mul(diag, tmp);
            free_mat(tmp);
            korrektes_urbild_ID(kernel_gen[i], nullcoz, G);
         }
         free_mat(nullcoz);

	 /*
         if (coz_info->aff_transl == NULL){
            coz_info->aff_transl = transl_aff_normal(G->gen, G->gen_no, &coz_info->aff_transl_no);
         }
         b1_kernel = kernel_factor_fct(coz_info->aff_transl, coz_info->aff_transl_no, G->gen_no, L,
                                       &b1_kernel_no);
         */

         /* calculate all preimages */
         all_preimages = calculate_all_preimages(preimage, kernel_gen, k, kernel_mat->cols,
                                                 H_GL_Z[1], b1_kernel, b1_kernel_no, anz);
         /* calculate the groups */
         SG = (bravais_TYP **)calloc(anz[0], sizeof(bravais_TYP *));
         for (i = 0; i < anz[0]; i++){
            SG[i] = extract_r(G, all_preimages[i]);
            free_mat(all_preimages[i]);
         }

         /* clean */
         free(all_preimages);
         free_mat(preimage);
	 /*
         for (i = 0; i < b1_kernel_no; i++){
            free_mat(b1_kernel[i]);
         }
         if (b1_kernel)
            free(b1_kernel);
	 */
         for (i = 0; i < kernel_mat->cols; i++){
            free_mat(kernel_gen[i]);
            free_mat(k[i]);
         }
         free(k);
         free(kernel_gen);
      }
      else{
         /* Berechne nur Vertreter unter dem affinen Normalisator */
	 /* ===================================================== */
	 reps = calculate_representatives(G, GL, H_G_Z, H_GL_Z, L, phi,
                                          H1_mod_ker, koord, kernel_mat, anz, length);
         SG = (bravais_TYP **)calloc(anz[0], sizeof(bravais_TYP *));
         for (i = 0; i < anz[0]; i++){
            tmp = darst_to_C1(H_GL_Z, reps[i]);
	    free_mat(reps[i]);
            tmp2 = mat_mul(diag, tmp);
            free_mat(tmp);
            korrektes_urbild(&tmp2, coz_info, G, L, FALSE);
	    SG[i] = extract_r(G, tmp2);
	    free_mat(tmp2);
	    length[0][i] *= b1_kernel_no;
	 }
	 free(reps);
      }
      for (i = 0; i < b1_kernel_no; i++){
         free_mat(b1_kernel[i]);
      }
      if (b1_kernel)
         free(b1_kernel);
   }
   else{
      SG = NULL;
   }

   /* clean */
   if (koord != NULL)
      free_mat(koord);
   free_mat(GLS);
   free_mat(diag);
   free_mat(kernel_mat);
   free_mat(phi);
   free_H1_mod_ker_TYP(H1_mod_ker);
   for (i = 0; i < image_gen_no; i++){
      free_mat(image[i]);
   }
   free(image);

   return(SG);
}







/* --------------------------------------------------------------------- */
/* Returns the maximal klassengleich subgroups of R with p_i-power-index */
/* where p_i is a prime given in G->divisors                             */
/* R has to be in standard affine form without translations              */
/* an G has to be the point group of R                                   */
/* with correct normalizer                                               */
/* Returns the number of these subgroups via anz.                        */

/* G->normal and G->gen have to generate the normalizer of G!!!!!        */

/* aflag: FALSE => calculate all max. k-subgroups                        */
/*        TRUE  => calculate max. k-subgroups up to conjugation of       */
/*                 the affine normalizer of R                            */
/* orbitlength: speichere die Laengen der Bahnen (nur wenn aflag)        */
/* --------------------------------------------------------------------- */
bravais_TYP **max_k_sub(bravais_TYP *G,
                        bravais_TYP *R,
                        matrix_TYP *pres,
                        int *anz,
			boolean aflag,
                        boolean debugflag,
			int **orbitlength)
{
   sublattice_TYP SL;

   bravais_TYP **S, *GL, **SG_L;

   matrix_TYP **H_G_Z, **H_GL_Z,
               *coz, *tmp, *diag, *preimage,
               *L, *Linv, *conj,
              **b1_kernel;

   word *relator;

   int i, j, k, sg_L_no, b1_kernel_no,
       size = MYSIZE, counter, *length = NULL;

   coz_TYP coz_info;

   rational eins;



   /* prepare */
   cen_to_norm(G);
   relator = (word *) calloc(pres->rows, sizeof(word));
   for (i = 0; i < pres->rows; i++){
      matrix_2_word(pres, relator + i, i);
   }
   coz = extract_c(R);
   eins.z = eins.n = 1;
   anz[0] = 0;
   S = (bravais_TYP **)calloc(0, sizeof(bravais_TYP *));
   if (aflag)
      orbitlength[0] = (int *)calloc(0, sizeof(int));

   /* calculate the sublattices */
   SL.lattices = max_sublattices(G, &SL.no, &SL.trivialflag, debugflag);

   /* calculate H^1(G, Q^n/Z^n) */
   H_G_Z = calculate_H1(G, relator, pres->rows);

   /* identify the cocyle */
   coz_info = identify_coz(G, coz, H_G_Z);

   /* calculate the orbits on the lattices under the stabilizer of the cocycle */
   SL.orbitlist = (int *)calloc(SL.no, sizeof(int));
   SL.orbitlength = (int *)calloc(SL.no + 1, sizeof(int));
   SL.smallest = (int *)calloc(SL.no + 1, sizeof(int));
   SL.conj = (matrix_TYP **)calloc(SL.no, sizeof(matrix_TYP *));
   SL.orbit_no = orbit_on_lattices(SL.lattices, SL.no, coz_info.Stab, coz_info.Stab_no,
                                   SL.orbitlist, SL.orbitlength, SL.smallest, SL.conj);

   for (i = 1; i <= SL.orbit_no; i++){
      /* calculate the maximal klassengleich subgroups for a representative of each orbit */
      L = SL.lattices[SL.smallest[i]];

      if (   TRUE   ){
/*
      if (!SL.trivialflag[SL.smallest[i]]){
*/
         /* lattice != p_i * Z^n */
         Linv = mat_inv(L);
         GL = my_konj_bravais(G, Linv, L);

         /* calculate H^1(GL, Q^n/Z^n) */
         H_GL_Z = calculate_H1(GL, relator, pres->rows);

         /* subgroups */
         SG_L = subgroups_for_L(&coz_info, L, G, GL, H_G_Z, H_GL_Z, &sg_L_no, &length, aflag);

         if (sg_L_no > 0){
            if (aflag){
               S = (bravais_TYP **)realloc(S, (anz[0] + sg_L_no)
                                            * sizeof(bravais_TYP *));
               orbitlength[0] = (int *)realloc(orbitlength[0], (anz[0] + sg_L_no)
                                            * sizeof(int));
            }
	    else{
               S = (bravais_TYP **)realloc(S, (anz[0] + sg_L_no * SL.orbitlength[i])
                                            * sizeof(bravais_TYP *));
            }
            for (j = 0; j < sg_L_no; j++){
               plus_translationen(SG_L[j], L);
               S[anz[0] + j] = SG_L[j];
	       if (aflag){
		  orbitlength[0][anz[0] + j] = length[j] * SL.orbitlength[i];
	       }
            }

            /* corresponding subgroups for the other lattices in this orbit */
	    if (!aflag){
               counter = 1;
               for (j = SL.smallest[i] + 1; j < SL.no; j++){
                  if (SL.orbitlist[j] == i){
                     conj = to_aff_normal_element(SL.conj[j], coz, 0, G, R);
                     for (k = 0; k < sg_L_no; k++){
                        S[anz[0] + counter * sg_L_no + k] =
                           konj_bravais(S[anz[0] + k], conj);
                     }
                     free_mat(conj);
                     counter++;
                  }
               }
               if (counter != SL.orbitlength[i]){ /* paranoiatest */
                  fprintf(stderr, "ERROR in max_k_sub!\n");
                  exit(8);
               }
               anz[0] += (sg_L_no * SL.orbitlength[i]);
	    }
	    else{
	       anz[0] += sg_L_no;
	    }
         }

         /* clean */
         for (j = 0; j < 3; j++)
            free_mat(H_GL_Z[j]);
         free(H_GL_Z);
         free_bravais(GL);
         free_mat(Linv);
         if (SG_L != NULL)
            free(SG_L);
         if (length)
	    free(length);
         length = NULL;
      }
      else{
         /* lattice == p_i * Z^n (there is only one lattice in this orbit) */
         /* calculate Ker(phi restricted to B^1(G,Q^n/L)) (here: phi = Id)*/

         /* this isn't correct!!! there has to be an additional condition
            example: min.56.1.1.1 min.56.1.1 2 3 */

	 /* Das ist nicht korrekt, weil | phi^{-1} (0 + B^1) |
	    nicht zwangslaeufig 1 ist. Dann wuerde es funktionieren! */

         if (coz_info.aff_transl == NULL){ /* we calculate it only once */
            coz_info.aff_transl = transl_aff_normal(G->gen, G->gen_no,
                                                    &coz_info.aff_transl_no);
         }
         b1_kernel = kernel_factor_fct(coz_info.aff_transl, coz_info.aff_transl_no,
                                       G->gen_no, L, &b1_kernel_no);
         if (size < (anz[0] + b1_kernel_no)){
            S = (bravais_TYP **)realloc(S, (anz[0] + b1_kernel_no) * sizeof(bravais_TYP *));
            size = anz[0] + b1_kernel_no;
         }

         if (coz_info.darst->cols == 0 || equal_zero(coz_info.darst)){
            preimage = init_mat(coz->rows, 1, "");
         }
         else{
            tmp = darst_to_C1(H_G_Z, coz_info.darst);
            diag = matrix_on_diagonal(L, G->gen_no);
            preimage = mat_mul(diag, tmp);
            free_mat(tmp);
            free_mat(diag);
         }
         korrektes_urbild(&preimage, &coz_info, G, L, TRUE);

         for (j = 0; j < b1_kernel_no; j++){
            b1_kernel[j] = mat_addeq(b1_kernel[j], preimage, eins, eins);
            S[j + anz[0]] = extract_r(G, b1_kernel[j]);
            free_mat(b1_kernel[j]);
            plus_translationen(S[j + anz[0]], L);
         }
         free(b1_kernel);
         free_mat(preimage);
         anz[0] += b1_kernel_no;
      }
      L = NULL;
   }

   /* paranoia test */
   if (debugflag){
      for (j = 0; j < anz[0]; j++){
         if (S[j] != NULL && !is_k_subgroup(S[j], R, G, G->gen[0]->rows, pres)){
            fprintf(stderr, "ERROR: das ist gar keine Untergruppe!\n");
            put_bravais(S[j],0,0);
            put_bravais(R,0,0);
            exit(8);
         }
      }
   }

   /* clean */
   free_sublattice_TYP(SL);
   for (i = 0; i < pres->rows; i++)
      wordfree(relator + i);
   free(relator);
   for (i = 0; i < 3; i++)
      free_mat(H_G_Z[i]);
   free(H_G_Z);
   free_coz_TYP(coz_info);
   free_mat(coz);

   return(S);
}




