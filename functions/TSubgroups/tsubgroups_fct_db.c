#include <typedef.h>
#include <base.h>
#include <bravais.h>
#include <matrix.h>
#include <presentation.h>
#include <tsubgroups.h>
#include <name.h>
#include <datei.h>


/* ---------------------------------------------------------------------- */
/* Berechne die maximalen t-Untergruppen einer Raumgruppe                 */
/* bis auf Konjugation unter der Raumgruppe unter Zuhilfename der         */
/* Datenbank.                                                             */
/* ---------------------------------------------------------------------- */
/* R: Raumgruppe in Standard affiner Form ohne Translationen              */
/* aflag: berechne Untergruppen bis auf Konjugation unter dem affinen     */
/*        Normalisator                                                    */
/* anzahl: speichere die Anzahl der Untergruppen hier                     */
/* ---------------------------------------------------------------------- */
TSubgroup_TYP **tsubgroup_db(bravais_TYP *R,
                             boolean aflag,
                             int *anzahl)
{
   CARATname_TYP Name;

   TSubgroup_TYP **sbg;

   bravais_TYP *Rstd, *Rinv;

   matrix_TYP **mat, *tmp, *inv;

   database *database;

   char pfad[1024], dbname[1024], format[1024];

   int i, k, dim, found = 0, j;



   /* berechne den Namen */
   dim = R->dim - 1;
   get_data_dir(dbname, "/tables/qcatalog/data");
   database = load_database(dbname, dim);
   Name = name_fct(R, database);


   /* suche Pfad (wird schon in name_fct bzw. deren Unterfunktionen
                  im Prinzip berechnet, koennte man da abspeichern) */
   i = 0;
   while (!found && i<database->nr ){
      found = (strcmp(database->entry[i].abbreviation, Name.qname) == 0);
      i++;
   }

   get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/");
   if (found){
      i--;
      sprintf(pfad, format,
              dim, database->entry[i].symbol,
              database->entry[i].order, database->entry[i].discriminant);

      Rstd = get_std_rep(pfad, Name);

      mat = get_words(pfad, Name, database->entry[i].affine, aflag, anzahl);

      sbg = (TSubgroup_TYP **)calloc(anzahl[0], sizeof(TSubgroup_TYP *));

      /* Nur die triviale Raumgruppe ist maximale Untergruppe! */
      if (anzahl[0] == 1 && mat[0]->rows == 1){
         sbg = (TSubgroup_TYP **)calloc(anzahl[0], sizeof(TSubgroup_TYP *));
         sbg[0] = (TSubgroup_TYP *)calloc(1, sizeof(TSubgroup_TYP));

         sbg[0]->R = init_bravais(R->dim);
         sbg[0]->R->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
         sbg[0]->R->gen_no = 1;
         sbg[0]->R->gen[0] = init_mat(R->dim, R->dim, "1");
         sbg[0]->orbitlength = 1;
         sbg[0]->pointgrouporder = 1;

         sbg[0]->P = init_bravais(R->dim - 1);
         sbg[0]->P->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
         sbg[0]->P->gen_no = 1;
         sbg[0]->P->gen[0] = init_mat(R->dim - 1, R->dim - 1, "1");
      }
      else{
         inv = mat_inv(Name.trafo);
         Rinv = init_bravais(Rstd->dim);
         Rinv->gen_no = Rstd->gen_no;
         Rinv->gen = (matrix_TYP **)calloc(Rstd->gen_no, sizeof(bravais_TYP *));
         for (k = 0; k < Rstd->gen_no; k++){
            Rinv->gen[k] = mat_inv(Rstd->gen[k]);
         }
         for (k = 0; k < anzahl[0]; k++){
            sbg[k] = ite_gruppe(Rstd, NULL, Rinv, NULL, mat[k], FALSE);
	    for (j = 0; j < sbg[k]->R->gen_no; j++){
               tmp = mat_kon(Name.trafo, sbg[k]->R->gen[j], inv);
	       free_mat(sbg[k]->R->gen[j]);
	       sbg[k]->R->gen[j] = tmp;
	       tmp = NULL;
	    }
         }
	 free_mat(inv);
	 free_bravais(Rinv);
      }
      
      if (mat){
         for (k = 0; k < anzahl[0]; k++)
            free_mat(mat[k]);
	 free(mat);
      }

      free_bravais(Rstd);
   }
   else{
      fprintf(stderr, "ERROR in tsubgroup_db!\n");
      exit(77);
   }

   /* clean */
   free_database (database);
   free_CARATname_TYP(Name);

   return(sbg);
}













