#include <typedef.h>
#include <base.h>
#include <bravais.h>
#include <matrix.h>
#include <presentation.h>
#include <tsubgroups.h>
#include <datei.h>
#include <graph.h>


/* -------------------------------------------------------------------- */
/* transform the i-th row of a matrix to a list of integers             */
/* -------------------------------------------------------------------- */
static int *row_to_list(matrix_TYP *mat,
                        int i)
{
   int *list, n;

   list = (int *)calloc(mat->cols, sizeof(int));
   for (n = 0; n < mat->cols; n++)
      list[n] = mat->array.SZ[i][n];

   return(list);
}

/* ---------------------------------------------------------------------- */
/* Setze R und P (nur falls aflag = TRUE), in die Zeilen von mat ein.     */
/* ---------------------------------------------------------------------- */
/* Rinv = R^{-1}                                                          */
/* Pinv = P^{-1}                                                          */
/* mat: Matrix mit Worten aus GAP                                         */
/*      (letzte Zeile: Laenge der Bahn unter Raumgr. und Pkt.gr.ordnung)  */
/* ---------------------------------------------------------------------- */
TSubgroup_TYP *ite_gruppe(bravais_TYP *R,
                          bravais_TYP *P,
			  bravais_TYP *Rinv,
			  bravais_TYP *Pinv,
                          matrix_TYP *mat,
			  boolean aflag)
{
   bravais_TYP *sg, *psg = NULL;

   int j, *word;

   TSubgroup_TYP *S;


   sg = init_bravais(R->dim);
   sg->gen_no = mat->rows - 1;
   sg->gen = (matrix_TYP **)calloc(sg->gen_no, sizeof(matrix_TYP *));
   if (aflag){
      psg = init_bravais(P->dim);
      psg->gen_no = mat->rows - 1;
      psg->gen = (matrix_TYP **)calloc(psg->gen_no, sizeof(matrix_TYP *));
   }

   for (j = 0; j < sg->gen_no; j++){
      word = row_to_list(mat, j);
      sg->gen[j] = mapped_word(word, R->gen, Rinv->gen);
      if (aflag)
         psg->gen[j] = mapped_word(word, P->gen, Pinv->gen);
      free(word);
   }

   /* in Struktur eintragen */
   /* ===================== */
   S = (TSubgroup_TYP *)calloc(1, sizeof(TSubgroup_TYP));
   S->R = sg;
   S->P = psg;
   S->orbitlength = mat->array.SZ[mat->rows - 1][0];
   S->pointgrouporder = mat->array.SZ[mat->rows - 1][1];
   return S;
}







/* ---------------------------------------------------------------------- */
/* Berechne die maximalen klassengleichen Untergruppen einer Raumgruppe   */
/* bis auf Konjugation unter der Raumgruppe.                              */
/* ---------------------------------------------------------------------- */
/* R: Raumgruppe in Standard affiner Form ohne Translationen              */
/* P: Punktgruppe von R (korrespondierende Erzeuger an der gleichen       */
/*    Position), P->gen, P->normal und P->cen zusammen muessen            */
/*    N_Gl_n(Z) (P) erzeugen.                                             */
/*    Der Formenraum muss korrekt sein.                                   */
/* pres: Praesentation von P                                              */
/* gapmats: Worte fuer die Konjugiertenklassen von maximalen Unter-       */
/*           gruppen von P aus GAP                                        */
/* no: Anzahl der Konjugiertenklassen                                     */
/* aflag: Raumgruppen bis auf Konjugation unter dem affinen Normalisator  */
/*        der Raumgruppe R                                                */
/* cflag: Ausgabe fuer Datenbank                                          */
/*        aflag wird auf true gesetzt                                     */
/* ---------------------------------------------------------------------- */
TSubgroup_TYP **tsubgroup(bravais_TYP *R,
                          bravais_TYP *P,
                          matrix_TYP *pres,
                          matrix_TYP **gapmats,
                          int *no,
                          boolean aflag,
			  boolean cflag)
{
   bravais_TYP *Rinv, *Pinv;

   TSubgroup_TYP **sbg, **S;

   int i, j, counter, *orbit,
       Nanz, Nanzalt, *factor;

   matrix_TYP **base, *M, **N;

   bahn **strong;

   CARATname_TYP Name;

   database *database;

   char dbname[1024];

   if (cflag)
      aflag = TRUE;
      
   /* lade Datenbank */
   get_data_dir(dbname, "/tables/qcatalog/data");
   database = load_database(dbname, P->dim);




   /* Nur die triviale Raumgruppe ist maximale Untergruppe! */
   /* ===================================================== */
   if (no[0] == 1 && gapmats[0]->rows == 1){
      sbg = (TSubgroup_TYP **)calloc(no[0], sizeof(TSubgroup_TYP *));
      sbg[0] = (TSubgroup_TYP *)calloc(1, sizeof(TSubgroup_TYP));

      sbg[0]->R = init_bravais(R->dim);
      sbg[0]->R->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      sbg[0]->R->gen_no = 1;
      sbg[0]->R->gen[0] = init_mat(R->dim, R->dim, "1");

      sbg[0]->orbitlength = 1;
      sbg[0]->pointgrouporder = 1;

      sbg[0]->P = init_bravais(P->dim);
      sbg[0]->P->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      sbg[0]->P->gen_no = 1;
      sbg[0]->P->gen[0] = init_mat(P->dim, P->dim, "1");

      if (cflag){
         Name = name_fct(R, database);
         printf("#1 %i %i ", Name.zname[0], Name.zname[1]);
         mpz_out_str(stdout, 10, &Name.aff_name);
	 printf("\n");
	 free_CARATname_TYP(Name);
         Name = name_fct(sbg[0]->R, database);
	 printf("0 1 %s %i %i ", Name.qname, Name.zname[0], Name.zname[1]);
         mpz_out_str(stdout, 10, &Name.aff_name);
	 printf("\n");
	 free_CARATname_TYP(Name);
      }
      free_database (database);
      return(sbg);
   }

   /* Die maximalen Untergruppen sind nicht trivial! */
   /* ============================================== */
   S = (TSubgroup_TYP **)calloc(no[0], sizeof(TSubgroup_TYP *));
   Rinv = init_bravais(R->dim);
   Rinv->gen_no = R->gen_no;
   Rinv->gen = (matrix_TYP **)calloc(R->gen_no, sizeof(bravais_TYP *));
   for (i = 0; i < R->gen_no; i++){
      Rinv->gen[i] = mat_inv(R->gen[i]);
   }
   if (aflag){
      Pinv = init_bravais(P->dim);
      Pinv->gen_no = P->gen_no;
      Pinv->gen = (matrix_TYP **)calloc(P->gen_no, sizeof(bravais_TYP *));
      for (i = 0; i < P->gen_no; i++){
         Pinv->gen[i] = mat_inv(P->gen[i]);
      }
   }

   for (i = 0; i < no[0]; i++){
      S[i] = ite_gruppe(R, P, Rinv, Pinv, gapmats[i], aflag);
   }

   free_bravais(Rinv);

   if (!aflag){
      free_database(database);
      return(S);
   }
   else{
      free_bravais(Pinv);

      /* berechne Punktgruppe des affinen Normalisators */
      /* ============================================== */
      cen_to_norm(P);
      N = PoaN(R, P, pres, &Nanz);
      Nanzalt = Nanz;
      Nanz = Nanz + P->gen_no;
      N = (matrix_TYP **)realloc(N, Nanz * sizeof(matrix_TYP *));
      for (i = 0; i < P->gen_no; i++)
         N[Nanzalt + i] = copy_mat(P->gen[i]);

      /* Fasse Bahnen unter R zu Bahnen unter affinem Normalisator zusammen */
      /* ================================================================== */
      orbit = (int *)calloc(no[0], sizeof(int));
      factor = (int *)calloc(no[0], sizeof(int));
      for (i = 0; i < no[0]; i++){
         orbit[i] = -1;
	 factor[i] = 0;
      }
      counter = 0;

      for (i = 0; i < no[0]; i++){
         if (orbit[i] == -1){
            orbit[i] = counter;
	    factor[counter] = 1;

            /* Strong generating set berechnen! */
            base = get_base(S[i]->P);
            strong = strong_generators(base, S[i]->P, TRUE);
            for (j = 0; j < P->dim; j++){
               free_mat(base[j]);
            }
            free(base);

            for (j = i + 1; j < no[0]; j++){
               if(S[i]->pointgrouporder == S[j]->pointgrouporder){

	          /* Die Funktion conjugated berechnet Bahn, so
		     dass eigentlich mehrere Gruppen auf einmal getestet
		     werden koennten. */
                  M = conjugated(S[i]->P, S[j]->P, N, Nanz, strong);
                  if (M != NULL){
                     orbit[j] = counter;
		     factor[counter]++;
                     free_mat(M);
                  }
               }
            }
            counter++;
            for (j = 0; j < P->dim; j++){
               free_bahn(strong[j]);
               free(strong[j]);
            }
            free(strong);
         }
      }

      /* Trage in Struktur ein */
      /* ===================== */
      sbg = (TSubgroup_TYP **)calloc(counter, sizeof(TSubgroup_TYP *));
      j = 0;
      if (cflag){
         Name = name_fct(R, database);
         printf("#%i %i %i ", counter, Name.zname[0], Name.zname[1]);
         mpz_out_str(stdout, 10, &Name.aff_name);
	 printf("\n");
	 free_CARATname_TYP(Name);
      }
      for (i = 0; i < no[0] && j < counter; i++){
         if (orbit[i] == j){
            sbg[j] = S[i];
	    sbg[j]->orbitlength *= factor[j];

	    if (cflag){
               Name = name_fct(sbg[j]->R, database);
	       printf("%i %i %s %i %i ", i, sbg[j]->orbitlength,
	                                 Name.qname, Name.zname[0], Name.zname[1]);
               mpz_out_str(stdout, 10, &Name.aff_name);
	       printf("\n");
	       free_CARATname_TYP(Name);
	    }

            j++;
	 }
         else{
            free_TSubgroup_TYP(S[i]);
         }
      }

      no[0] = counter;
      free(orbit);
      free(factor);
      free_database(database);
      free(S);
      for (i = 0; i < Nanz; i++)
         free_mat(N[i]);
      free(N);
      return(sbg);
   }
}















