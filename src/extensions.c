#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <datei.h>
#include <gmp.h>
#include <zass.h>
#include <longtools.h>
#include <presentation.h>
#include <base.h>
#include <graph.h>

#define DEBUG FALSE

int INFO_LEVEL;
extern int SFLAG;


static int is_trivial(bravais_TYP *G){

  int i;

  for (i=0;i<G->gen_no;Check_mat(G->gen[i]),i++);
  for (i=0;i<G->gen_no &&
           G->gen[i]->flags.Scalar &&
           G->gen[i]->array.SZ[i][i] == 1;i++);

  return i == G->gen_no;

}

int main(int argc,char **argv){

  matrix_TYP **X,
              *N,
             **Y,
             **TR,
             **base,
              *relator_input,
             **matinv;

  word *relator;

  bravais_TYP *G, *R;

  int i,
      anz,
     *len;

  long dim;

  bahn **strong;


  MP_INT number,
        *names;

  char comment[1000],
       file[1000],
      *NAME,
      *FN;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)
      || (is_option('i') && FILEANZ<2)){
      printf("Usage: %s ['file1'] 'file2' ['file3'] [-n] [-i] [-t=n] [-v] [-C] [-H] [-F] [-S]\n",argv[0]);
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
      printf(" -t=n:    Has an effect only if given with -i. Outputs the\n");
      printf("          isomophism needed to transform the space group\n");
      printf("          By default, only the linear part is calculated. To\n");
      printf("          get a full transformation matrix, use -t=2.\n");
      printf(" -C:      Ignore the operation of the normalizer, just work\n");
      printf("          on the level of extensions.\n");
      printf(" -H:      echo the isomorphism type of the cohomology group\n");
      printf("          H^1(G,Q^n/Z^n) to stderr.\n");
      printf(" -F:      Only construct those extensions which gie rise to\n");
      printf("          torsion free space groups. Does not work in conjunction\n");
      printf("          with -n.\n");
      printf(" -S:      write corresponding space groups in files\n");
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


  if ((FILEANZ > 2)
      || (FILEANZ == 2 && ! is_option('i'))){
    /* reading the group */
    G = get_bravais(FILENAMES[1]);

    /* reading the presentation */
    X = mget_mat(FILENAMES[0],&i);
    if (i>1){
      fprintf(stderr,"you should only give a single matrix as presentation\n");
      exit(3);
    }

    FN = FILENAMES[1];
  }
  else{
    /* reading the group */
    G = get_bravais(FILENAMES[0]);

    /* calculate a presentation */
    base = get_base(G);

    strong = strong_generators(base,G,TRUE);

    if (DEBUG){
      check_base(strong,G);
    }

    X = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
    X[0] = pres(strong,G,NULL);

    for (i=0;i<G->dim;i++){
      free_mat(base[i]);
      free_bahn(strong[i]);
      free(strong[i]);
    }
    free(strong);
    free(base);
    FN = FILENAMES[0];
  }


  relator_input = X[0];
  free(X);

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
        printf("number of extensions of group in %s with natural lattice %d\n",
              FN,1);
     }
     else if(is_option('i')){
        printf("There is only one extension of this group with the natural\n");
        printf("lattice, and this splits.\n");
     }
     else if(is_option('F')){
        if (is_trivial(G)){
           X[0] = init_mat(G->gen_no * G->dim,1,"");
           printf("#%d\n",1);
           sprintf(comment,"%% the %d-th cozycle to the group of %s",
                   1,FN);
           if (is_option('S')){
              sprintf(file, "%s.%s", FN, "0");
              R = extract_r(G, X[0]);
	      put_bravais(R, file, comment);
	      free_bravais(R);
	   }
	   else{
              put_cocycle(X[0],G->dim,G->gen_no,NULL,comment);
	   }
           free_mat(X[0]);
        }
        else{
           /* only one extension, this splits and is not torsion free
              for this reason */
           printf("#0\n");
        }
     }
     else{
        X[0] = init_mat(G->gen_no * G->dim,1,"");
        printf("#%d\n",1);
        sprintf(comment,"%% the %d-th cozycle to the group of %s",
                1,FN);
        if (is_option('S')){
           sprintf(file, "%s.%s", FN, "0");
           R = extract_r(G, X[0]);
	   put_bravais(R, file, comment);
	   free_bravais(R);
	}
	else{
           put_cocycle(X[0],G->dim,G->gen_no,NULL,comment);
	}
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
              FN);
      mpz_out_str(stdout,10,&number);
      printf("\n");
      mpz_clear(&number);
  }
  else if(is_option('i')){

     if (FILEANZ > 2){
       Y = mget_mat(FILENAMES[2],&anz);
     }
     else{
       Y = mget_mat(FILENAMES[1],&anz);
     }

     convert_cocycle_to_column(Y,anz,G->dim,G->gen_no);
     names = (MP_INT *) malloc(anz*sizeof(MP_INT));
     for (i=0;i<anz;i++) mpz_init_set_si(names+i,0);
     i = is_option('t');
     if (i && optionnumber('t') == 2) i = 3;
     TR = identify(X[0],X[1],X[2],G,Y,names,anz,i,NULL,NULL);
     for (i=0;i<anz;i++){
        printf("Name for the %d-th extension in %s: ",i+1,FN);
        mpz_out_str(stdout,10,names+i);
        printf("\n");
        if (TR!=NULL){
           sprintf(comment,"transformation matrix to extension %d of %s",
                    i+1,FN);
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
     Y = extensions(X[0],X[1],X[2],G,&len,&names,&anz,is_option('F'));

     printf("#%d\n",anz);
     for (i=0;i<anz;i++){
        NAME = mpz_get_str(NULL,10,names+i);
        sprintf(comment,
             "%% the %d-th cozycle, length of orbit %d,name: %s",i+1,len[i],NAME);
        if (is_option('S')){
           sprintf(file, "%s.%s", FN, NAME);
           R = extract_r(G, Y[i]);
	   put_bravais(R, file, comment);
	   free_bravais(R);
	}
	else{
           put_cocycle(Y[i],G->dim,G->gen_no,NULL,comment);
	}
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

