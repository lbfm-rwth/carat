#include <typedef.h>
#include <bravais.h>
#include <name.h>
#include <tsubgroups.h>
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include <string.h>
#include "datei.h"


/* -------------------------------------------------------- */
/* Free TSubgroup_TYP and pointer to this!                  */
/* -------------------------------------------------------- */
void free_TSubgroup_TYP(TSubgroup_TYP *sbg)
{
   if (sbg){
      if (sbg->R)
         free_bravais(sbg->R);
      if (sbg->P)
         free_bravais(sbg->P);
      free(sbg);
   }
}



/* -------------------------------------------------------- */
/* Free CARATname_TYP!                                      */
/* -------------------------------------------------------- */
void free_CARATname_TYP(CARATname_TYP Name)
{
   if (&Name.aff_name)
      mpz_clear(&Name.aff_name);
   if (Name.trafo)
      free_mat(Name.trafo);
}




/* -------------------------------------------------------- */
/* Free TSUB_TYP!                                           */
/* -------------------------------------------------------- */
void free_TSUB_TYP(TSUB_TYP TSUB)
{
   int i;
   
   for (i = 0; i < TSUB.word_no; i++){
      if (TSUB.words[i])
         free_mat(TSUB.words[i]);
   }
}


/* -------------------------------------------------------- */
/* Get standardrepresentative of an affine class            */
/* (compare with reverse_name_fct.c)                        */
/* -------------------------------------------------------- */
/* pfad: Pfad, wo Q-Klassen-Vertreter ist                   */
/* Name: CARAT-Name                                         */
/* -------------------------------------------------------- */
bravais_TYP *get_std_rep(char *pfad,
                         CARATname_TYP Name)
{
   bravais_TYP *R;

   char filename[2048];

   bravais_TYP *DATAZ, *DATAQ;

   matrix_TYP *PRES;


   /* calculate standard representative */
   sprintf(filename, "%s/%s", pfad, Name.qname);
   DATAQ = get_bravais(filename);
   sprintf(filename,"%s/pres.%s", pfad, Name.qname);
   PRES = get_mat(filename);
   DATAZ = get_zclass_by_name(DATAQ, Name.zname, Name.zname+1, FALSE);
   R = get_affine_class_by_name(DATAZ, PRES, &Name.aff_name, 1);

   /* clean */
   cleanup_prime();
   free_bravais(DATAQ);
   free_bravais(DATAZ);
   free_mat(PRES);

   return(R);
}




/* -------------------------------------------------------- */
/* hole die Worte fuer die t-Untergruppen einer Raumgruppe  */
/* aus der Datenbank                                        */
/* -------------------------------------------------------- */
/* pfad: Verzeichnis der Q-Klasse zu der Raumgruppe         */
/* Name: CARATname der Raumgruppe                           */
/* aff_class_no: Anzahl der aff. Klassen in der Q-Klasse    */
/* aflag: berechne die t-Untergr. bis auf Konjugation unter */
/*        dem affinen Normalisator                          */
/* anzahl: speichere die Anzahl der t-Untergr. hier         */
/* -------------------------------------------------------- */
matrix_TYP **get_words(char *pfad,
                       CARATname_TYP Name,
		       int aff_class_no,
		       boolean aflag,
		       int *anzahl)
{
  int anz, k, i, z1, z2, *woerter, nr, laenge, c;

   char filename[2048], string[512];

   matrix_TYP **mat, **tmp;

   FILE *infile;

   MP_INT aff_name;

   boolean FLAG;


   /* oeffne Datei */
   sprintf(filename, "%s/words.%s", pfad, Name.qname);
   if ( (infile = fopen(filename, "r")) == NULL ) {
      fprintf(stderr, "get_words: Error: Could not open input-file!\n");
      exit (4);
   }

   /* Hole Worte fuer alle Untergruppen */
   c=fscanf (infile, "%[^\n]",string);
   if ( string[0] != '#' ) {
      anz = 1;
      mat = (matrix_TYP **)malloc(sizeof(matrix_TYP *));
      rewind(infile);
      mat[0] = fget_mat(infile);
   }
   else{
      sscanf (string, "#%u", &anz);
      mat = (matrix_TYP **)malloc(anz * sizeof(matrix_TYP *));
      for (k = 0; k < anz; k++){
         mat[k] = fget_mat(infile);
      }
   }

   /* trivialer Fall */
   if (anz == 0){
      anzahl[0] = 0;
      fclose(infile);
      return(NULL);
   }

   if (aflag){
      woerter = (int *)calloc(anz, sizeof(int));
      FLAG = FALSE;

      /* suche die richtige affine Klasse */
      for (k = 0; k < aff_class_no; k++){
         if (fscanf(infile, "%s%i%i", string, &z1, &z2) != 3){
            fprintf (stderr, "Data has wrong structure.\n");
            exit (4);
	 }
         if (sscanf (string, "#%i", anzahl) != 1){
            fprintf (stderr, "Data has wrong structure.\n");
            exit (4);
         }
         mpz_init(&aff_name);
         mpz_inp_str(&aff_name, infile, 10);

	 if (z1 == Name.zname[0] && z2 == Name.zname[1] &&
	     mpz_cmp(&aff_name, &Name.aff_name) == 0){
	    FLAG = TRUE;
	 }
	 mpz_clear(&aff_name);

	 /* lese Infos ueber die Vertreter der t-Untergr. bis auf Konj.
	    unter dem aff. Normalisator */
         for (i = 0; i < anzahl[0]; i++){
            if (fscanf(infile, "%i%i%s%i%i", &nr, &laenge, string, &z1, &z2) != 5){
               fprintf (stderr, "Data has wrong structure.\n");
               exit (4);
            }
            mpz_init(&aff_name);
            mpz_inp_str(&aff_name, infile, 10);
	    mpz_clear(&aff_name);
            if (FLAG){
	       woerter[nr]++;

	       /* Laenge der Bahn unter dem affinen Normalisator */
	       mat[nr]->array.SZ[mat[nr]->rows - 1][0] = laenge;
	    }
	 }
	 if (FLAG)
	    break;
      }
      if (!FLAG){
         fprintf(stderr, "ERROR in get_words!\n");
	 exit(56);
      }

      /* lasse nur die Matrizen fuer die Vertreter uebrig */
      tmp = (matrix_TYP **)calloc(anz, sizeof(matrix_TYP *));
      i = 0;
      for (k = 0; k < anz; k++){
         if (woerter[k]){
	    if (woerter[k] != 1){
	       fprintf(stderr, "ERROR in  get_words!\n");
	       exit(58);
	    }
	    tmp[i] = mat[k];
	    i++;
	 }
	 else
	    free_mat(mat[k]);
      }
      if (i != anzahl[0]){
         fprintf(stderr, "ERROR in get_words!\n");
	 exit(57);
      }
      free(mat);
      mat = tmp;
      tmp = NULL;
      free(woerter);
   }
   else{
      anzahl[0] = anz;
   }

   fclose(infile);

   return(mat);
}



















