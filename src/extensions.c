#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>

#include <gmp.h>
#include <zass.h>
#include <longtools.h>


int INFO_LEVEL;
extern int SFLAG;

static void put_cozycle(matrix_TYP *COZ,
                        int dim,
                        char *file,
                        char *comment)
{
   int i,
       rr = 0;

   FILE *F;

   if (file == NULL){
      F = stdout;
   }
   else{
      F = fopen(file,"rw");
   }

   if (F == NULL){
      fprintf(stderr,"problems opening %s\n",file);
      exit(4);
   }

   if (COZ->cols != 1){
      fprintf(stderr,"cozycle with more than 1 columns?\n");
      exit(3);
   }

   rat2kgv(COZ);
   Check_mat(COZ);

   if (COZ->kgv == 1 ||
       COZ->kgv == 0){
      fprintf(F,"%dx%d\t%s\n",COZ->rows,COZ->cols,comment);
   }
   else{
      fprintf(F,"%dx%d/%d\t%s\n",COZ->rows,COZ->cols,COZ->kgv,comment);
   }

   for (i=0;i<COZ->rows;i++){
      
      fprintf(F,"%d ",COZ->array.SZ[i][0]);

      rr = (rr+1) % dim;
      if (rr == 0){
         fprintf(F,"\n");
      }
   }

   return;
}



main(int argc,char **argv){

  matrix_TYP **X,
              *N,
             **Y,
             **TR,
              *relator_input,
             **matinv;

  word *relator;

  bravais_TYP *G;

  int i,
      anz,
     *len;

  long dim;

  MP_INT number,
        *names;

  char comment[1000],
      *NAME;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 2)
      || (is_option('i') && FILEANZ<3)){
      printf("Usage: %s 'file1' 'file2' ['file3'] [-n] [-i] [-t] [-v]\n",argv[0]);
      printf("\n");
      printf("file1:  matrix_TYP containing a presentation of the group (cf. Presentation,\n");
      printf("        Roundcor)\n");
      printf("file2:  bravais_TYP containing the group (and its normalizer).\n");
      printf("file3:  (OPTIONAL) containing cocycles for identification.\n");
      printf("\n");
      printf("Calculates representatives of the isomorphism classes of extension of\n");
      printf("the group in file2 with presentation in file1 with the natural lattice.\n");
      printf("Output lists representative vectorsystems. (Output can be used as input\n");
      printf("for Extract). Also the isomorphism type of the extension group is given.\n");
      printf("\n");
      printf("Options:\n");
      printf("\n");
      printf(" -n:      only calculates the number of extensions, without computing\n");
      printf("          representatives. This option is much faster for big cohomology\n");
      printf("          groups. WARNING: Try this option first for 2-groups in dimension\n");
      printf("          greater than 4.\n");
      printf(" -v:      verbose mode. Give some echoing to stderr to\n");
      printf("          indicate a little what the program is doing.\n");
      printf(" -i:      Identify the cozycles given in file3, ie. give\n");
      printf("          the described space groups a name. CAUTION: The\n");
      printf("          name will depend on the generating set of the group\n");
      printf("          in file2 and the presentation in file1.\n");
      printf("          Can be used to test isomorphism of space groups with equal\n");
      printf("          point groups. The name is 0 iff the extension splits.\n");
      printf(" -t:      Has an effect only if given with -i. Outputs the\n");
      printf("          isomophism needed to transform the space group\n");
      printf(" -C:      Ignore the operation of the normalizer, just work\n");
      printf("          on the level of extensions.\n");
      printf(" -H:      echo the isomorphism type of the cohomology group\n");
      printf("          H^1(G,Q^n/Z^n) to stderr.\n");
      printf("\n");
      printf(" CAUTION: THE PROGRAM RELIES HEAVILY ON THE FOLLOWING TWO FACTS:\n");
      printf("           - THE PRESENTATION GIVEN IN file1 IS INDEED A\n");
      printf("             PRESENTATION. CHECKED IS ONLY THAT THE GIVEN\n");
      printf("             GENERATORS FULFILL THE RELATIONS.\n");
      printf("           - THE NORMALIZER IN GL_n(Z) FOR THE GROUP IN file2\n");
      printf("             IS CORRECT. THE PROGRAM MIGHT NOT FINISH OR\n");
      printf("             GIVE WRONG ANSWERS OTHERWISE.\n");
      printf("\n");
      printf("Cf. also Zass_main, Extract, Same_generators, Presentation.\n");
      printf("Synonyms: Vectorsystems,Extensions\n");
      if (is_option('h')){
         exit(0);
      }
      else{
         exit(31);
      }
  }

  INFO_LEVEL = optionnumber('h');

  if (INFO_LEVEL & 12){
     SFLAG = 1;
  }

  /* reading the presentation */
  X = mget_mat(FILENAMES[0],&i);
  if (i>1){
    fprintf(stderr,"you should only give a single matrix as presentation\n");
    exit(3);
  }
  relator_input = X[0];
  free(X);

  /* reading the group */
  G = get_bravais(FILENAMES[1]);

  /* if C is an option, completely forget the normalizer and centralizer */
  if (is_option('C')){
     for (i=0;i<G->normal_no;i++) free_mat(G->normal[i]);
     if (G->normal_no > 0){
        G->normal_no = 0;
        free(G->normal);
     }
     for (i=0;i<G->cen_no;i++) free_mat(G->cen[i]);
     if (G->cen_no > 0){
        G->cen_no = 0;
        free(G->cen);
     }
  }

  /* we have to have at least the identity to generate the normalizer */
  if (G->normal_no == 0){
     G->normal_no = 1;
     G->normal = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP *));
     G->normal[0] = init_mat(G->dim,G->dim,"1");
  }

  /* bereitstellen eines pointers fuer die inversen von mat */
  matinv = (matrix_TYP **) calloc(G->gen_no,sizeof(matrix_TYP *));

  /* speicher fuer die worte */
  relator = (word *) calloc(relator_input->rows,sizeof(word));

  /* konvertieren der inputmatrix in relator-format */
  for (i=0;i<relator_input->rows;i++){
    matrix_2_word(relator_input,relator+i,i);
  }

  X = cohomology(&dim,G->gen,matinv,relator,G->gen_no,relator_input->rows);

  /* there is a special case to handle, which is the case that there
     isn't a cohomology group at all */
  if (X[0]->cols <1){
     if (is_option('H')){
        fprintf(stderr,"H^1(G,Q^n/Z^n) is trivial\n");
     }
     if (is_option('n')){
        printf("number of extensions of group in %s with natural lattice %d",
              FILENAMES[1],1);
     }
     else if(is_option('i')){
        printf("There is only one extension of this group with the natural\n");
        printf("lattice, and this splits.\n");
     }
     else{
        X[0] = init_mat(G->gen_no * G->dim,1,"");
        printf("#%d\n",1);
        sprintf(comment,"the %d-th cozycle to the group of %s",
                1,FILENAMES[1]);
        put_mat(X[0],NULL,comment,2);
     }
     exit(0);
  }

  if (is_option('H')){
     fprintf(stderr,"Isomorphism type of H^1(G,Q^n/Z^n): \n");
     i=0;
     while(X[1]->array.SZ[i][i] == 1) i++;
     fprintf(stderr,"C_%d",X[1]->array.SZ[i][i]);
     for (i=i+1;i<X[1]->cols;i++){
        if (X[1]->array.SZ[i][i] != 0 &&
            X[1]->array.SZ[i][i] != 1){
           fprintf(stderr," X ");
           fprintf(stderr,"C_%d",X[1]->array.SZ[i][i]);
        }
     }
     fprintf(stderr,"\n");
  }

  if (is_option('n')){
      mpz_init(&number);
      mpz_set_si(&number,0);
      no_of_extensions(X[0],X[1],X[2],G,&number);
      printf("number of extensions of group in %s with natural lattice ",
              FILENAMES[1]);
      mpz_out_str(stdout,10,&number);
      printf("\n");
      mpz_clear(&number);
  }
  else if(is_option('i')){
     Y = mget_mat(FILENAMES[2],&anz);
     names = (MP_INT *) malloc(anz*sizeof(MP_INT));
     for (i=0;i<anz;i++) mpz_init_set_si(names+i,0);
     i = is_option('t');
     TR = identify(X[0],X[1],X[2],G,Y,names,anz,i);
     for (i=0;i<anz;i++){
        printf("Name for the %d-th extension in %s: ",i+1,FILENAMES[1]);
        mpz_out_str(stdout,10,names+i);
        printf("\n");
        if (TR!=NULL){
           sprintf(comment,"transformation matrix to extension %d of %s",
                    i+1,FILENAMES[1]);
           put_mat(TR[i],NULL,comment,2);
        }
     }
     for (i=0;i<anz;i++){
        mpz_clear(names+i);
        free_mat(Y[i]);
        if (TR != NULL) free_mat(TR[i]);
     }
     free(names);
     free(Y);
     if (TR != NULL) free(TR);
  }
  else{
     Y = extensions(X[0],X[1],X[2],G,&len,&names,&anz);

     printf("#%d\n",anz);
     for (i=0;i<anz;i++){
        NAME = mpz_get_str(NULL,10,names+i);
        sprintf(comment,
             "the %d-th cozycle, length of orbit %d,name: %s",i+1,len[i],NAME);
        put_cozycle(Y[i],G->dim,NULL,comment);
        free_mat(Y[i]);
        mpz_clear(names+i);
        free(NAME);
     }
     free(names);
     free(Y);
     Y = NULL;
     free(len);
  }

  for (i=0;i<3;i++) free_mat(X[i]);
  for (i=0;i<G->gen_no;i++) if (matinv[i] != NULL) free_mat(matinv[i]);
  free(matinv);
  for (i=0;i<relator_input->rows;i++) wordfree(relator+i);
  free(relator);
  free_mat(relator_input);
  free(X);
  free_bravais(G);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }


  exit(0);
} /* main */

