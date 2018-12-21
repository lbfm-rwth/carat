#include <typedef.h>
#include <base.h>
#include <bravais.h>
#include <matrix.h>
#include <getput.h>
#include <presentation.h>
#include <tsubgroups.h>
#include <name.h>
#include <datei.h>


/* ---------------------------------------------------------------------- */
/* Berechne diejenigen t-Obergruppen einer Raumgruppe, die in einer       */
/* gegebenen Q-Klasse sind.                                               */
/* ---------------------------------------------------------------------- */
/* pfad: Pfad zur Q-Klasse, die ueberprueft werden soll                   */
/* qnameS: Name der Q-Klasse                                              */
/* RName: Name der affinen Klasse der Raumgruppe                          */
/* anzahl: speichere die Anzahl der Obergruppen hier                      */
/* aff_class_no: Anzahl der affinen Klassen in der Q-Klasse               */
/* database: Datenbank der Q-Klassen in einer Dimension                   */
/* trafoinv: RName.trafo^{-1}                                             */
/* ---------------------------------------------------------------------- */
static bravais_TYP **get_supergr(char *pfad,
                                 char *qnameS,
                                 CARATname_TYP RName,
                                 int *anzahl,
				 int aff_class_no,
				 database *database,
				 matrix_TYP *trafoinv)
{
  int i, k, j, z1ri, z2ri, number, anz, nr, laenge, c;

   bravais_TYP *Sstd, *Ri, **Si = NULL, *Sstdinv;

   CARATname_TYP NameSstd, NameRi;

   char filename[1024], string[512];

   matrix_TYP **mat, *inv, *trafo;

   FILE *infile;

   MP_INT aff_name_ri;

   TSubgroup_TYP *SG;




   /* Vorbereitungen */
   anzahl[0] = 0;
   sprintf(filename, "%s/words.%s", pfad, qnameS);
   if ( (infile = fopen(filename, "r")) == NULL ) {
      fprintf(stderr, "get_supergr: Error: Could not open input-file!\n");
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
      fclose(infile);
      return(NULL);
   }

   /* suche die affinen Klassen, die unsere affine Klasse als Untergruppe haben */
   for (k = 0; k < aff_class_no; k++){ /* durchlaufe alle affinen Klassen */
      Sstd = NULL;
      strcpy(NameSstd.qname, qnameS);
      if (fscanf(infile, "%s%i%i", string, &NameSstd.zname[0], &NameSstd.zname[1]) != 3){
         fprintf (stderr, "get_supergr: Error: Data has wrong structure.\n");
         exit (4);
      }
      if (sscanf (string, "#%i", &number) != 1){
         fprintf (stderr, "get_supergr: Error: Data has wrong structure.\n");
         exit (4);
      }
      mpz_init(&NameSstd.aff_name);
      mpz_inp_str(&NameSstd.aff_name, infile, 10);

      /* lese Infos ueber die Vertreter der t-Untergr. bis auf Konj.
	 unter dem aff. Normalisator */
      for (i = 0; i < number; i++){ /* durchlaufe alle Untergruppenbahnen der aff. Kl. */
         if (fscanf(infile, "%i%i%s%i%i", &nr, &laenge, string, &z1ri, &z2ri) != 5){
            fprintf (stderr, "get_supergr: Error: Data has wrong structure.\n");
            exit (4);
         }
         mpz_init(&aff_name_ri);
         mpz_inp_str(&aff_name_ri, infile, 10);

	 if (strcmp(RName.qname, string) == 0 &&
	     z1ri == RName.zname[0] && z2ri == RName.zname[1] &&
	     mpz_cmp(&aff_name_ri, &RName.aff_name) == 0){
            /* wir haben eine solche affine Klasse gefunden */

	    if (RName.order == 1){ /* die Raumgruppe ist die trivial */
	       if (anzahl[0] == 0)
                  Si = (bravais_TYP **)calloc(1, sizeof(bravais_TYP *));
	       else
                  Si = (bravais_TYP **)realloc(Si, (anzahl[0] + 1) * sizeof(bravais_TYP *));
               Si[anzahl[0]] = get_std_rep(pfad, NameSstd);
               anzahl[0]++;
	    }
	    else{ /* die Raumgruppe ist nicht trivial */

	       /* berechne Standardvertreter der Obergruppe */
               if (Sstd == NULL){
	          Sstd = get_std_rep(pfad, NameSstd);
                  Sstdinv = init_bravais(Sstd->dim);
                  Sstdinv->gen_no = Sstd->gen_no;
                  Sstdinv->gen = (matrix_TYP **)calloc(Sstd->gen_no, sizeof(bravais_TYP *));
                  for (j = 0; j < Sstd->gen_no; j++)
                     Sstdinv->gen[j] = mat_inv(Sstd->gen[j]);
               }

               /* berechne die Untergruppe von Sstd */
               SG = ite_gruppe(Sstd, NULL, Sstdinv, NULL, mat[nr], FALSE);
               Ri = SG->R;
	       SG->R = NULL;
	       free_TSubgroup_TYP(SG);

	       /* berechne Matrix, die Ri zu Standardvertreter Rstd konjugiert */
	       /* (dies koennte man wesentlich eleganter machen, da man ja den
                  Standardvertreter schon kennt und ihn nicht erst suchen muss) */
               NameRi = name_fct(Ri, database);
               trafo = mat_mul(NameRi.trafo, trafoinv);
	       inv = mat_inv(trafo);
	       free_CARATname_TYP(NameRi);

	       /* konjugiere Sstd zu Obergruppe von Standardvertreter Rstd */
	       if (anzahl[0] == 0)
                  Si = (bravais_TYP **)calloc(1, sizeof(bravais_TYP *));
	       else
                  Si = (bravais_TYP **)realloc(Si, (anzahl[0] + 1) * sizeof(bravais_TYP *));
	       Si[anzahl[0]] = init_bravais(Sstd->dim);
	       Si[anzahl[0]]->gen_no = Sstd->gen_no;
	       Si[anzahl[0]]->gen = (matrix_TYP **)calloc(Sstd->gen_no, sizeof(matrix_TYP *));
	       for (j = 0; j < Sstd->gen_no; j++){
	          Si[anzahl[0]]->gen[j] = mat_kon(inv, Sstd->gen[j], trafo);
	       }
               anzahl[0]++;
	       free_mat(inv);
	       free_mat(trafo);
	       free_bravais(Ri);
	    }
	 }
	 mpz_clear(&aff_name_ri);
      }
      if (Sstd){
         free_bravais(Sstd);
	 free_bravais(Sstdinv);
      }
      free_CARATname_TYP(NameSstd);
   }
   
   /* clean */
   for (k = 0; k < anz; k++)
      free_mat(mat[k]);
   free(mat);


   return(Si);
}




/* ---------------------------------------------------------------------- */
/* Berechne die minimalen translationengleichen Obergruppen einer         */
/* Raumgruppe bis auf Konjugation unter dem affinen Normalisator der      */
/* Raumgruppe unter Zuhilfename  der Datenbank.                           */
/* ---------------------------------------------------------------------- */
/* R: Raumgruppe in Standard affiner Form ohne Translationen              */
/* anzahl: speicher die Anzahl der Raumgruppen dort                       */
/* ---------------------------------------------------------------------- */
bravais_TYP **tsupergroups(bravais_TYP *R,
                           int *anzahl)
{
   CARATname_TYP Name;

   bravais_TYP **Super, **tmp;

   matrix_TYP *inv;

   database *database;

   char pfad[1024], dbname[1024], format[1024];

   int i, j, no, dim;


   /* Vorbereitungen */
   anzahl[0] = 0;
   Super = (bravais_TYP **)calloc(0, sizeof(bravais_TYP *));

   /* lade Datenbank */
   dim = R->dim - 1;
   get_data_dir(dbname, "/tables/qcatalog/data");
   database = load_database(dbname, dim);

   /* berechne den Namen */
   Name = name_fct(R, database);
   inv = mat_inv(Name.trafo);

   get_data_dir(format, "tables/qcatalog/dim%d/dir.%s/ordnung.%d/%s/");
   for (i = 0; i < database->nr; i++){
      if (database->entry[i].order > Name.order &&
          database->entry[i].order % Name.order == 0){
	  sprintf(pfad, format,
                 dim, database->entry[i].symbol,
                 database->entry[i].order, database->entry[i].discriminant);

         tmp = get_supergr(pfad, database->entry[i].abbreviation, Name,
			   &no, database->entry[i].affine, database, inv);

	 if (no > 0){
	    Super = (bravais_TYP **)realloc(Super, (anzahl[0] + no) * sizeof(bravais_TYP *));
  	    for (j = 0; j < no; j++){
               Super[anzahl[0] + j] = tmp[j];
	       tmp[j] = NULL;
	    }
	    free(tmp);
            anzahl[0] += no;
	 }
      }
   }


   /* clean */
   free_database (database);
   free_mat(inv);
   free_CARATname_TYP(Name);

   return(Super);
}













