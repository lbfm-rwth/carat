#include <typedef.h>

#include <getput.h>
#include <matrix.h>
#include <orbit.h>
#include <bravais.h>
#include <base.h>
#include <tools.h>
#include <zass.h>
#include <datei.h>
#include <presentation.h>

int INFO_LEVEL;
extern int SFLAG;

matrix_TYP *translation_lattice(matrix_TYP **G,int number,matrix_TYP *P);

int main(int argc,char **argv){

  bravais_TYP *G,
              *H;

  matrix_TYP *X,
            **XX,
            **base,
             *T,
            **COZ;

  bahn **strong;

  int i,
      j,
      k,
      type,         /* TRUE iff the function has been called via 
                        ...Standart_affine_form */
      kgv,
      cozanz;

  char comment[1000],
       file[1000];

  G = NULL;
  H = NULL;
  X = NULL;
  T = NULL;

  read_header(argc,argv);

  #define OTHER_NAME "Standard_affine_form"
  if (strlen(argv[0]) < strlen(OTHER_NAME)){
     /* it can't have been called via .....OTHER_NAME */
     type = FALSE;
  }
  else{
     type = (strcmp(argv[0]+(strlen(argv[0])-strlen(OTHER_NAME)),
             OTHER_NAME) == 0);
  }


  if ((type && FILEANZ < 1) ||
      (!type && ((is_option('h') && optionnumber('h')==0)||
      (FILEANZ < 1) ||
      ((FILEANZ < 2) && (is_option('r')))))){
     if (type){
        printf("Usage: %s 'file1' 'file2' [-t]\n",argv[0]);
        printf("\n");
        printf("file1: bravais_TYP containing a space group G.\n");
        printf("file2: matrix_TYP contatining a presentation of the point\n");
        printf("       group of G.\n");
        printf("\n");
        printf("Transforms the space group with generators in  file1 into standard\n");
        printf("form, ie. the translation lattice is transformed to Z^n.\n");
        printf("In case the translation subgroup has dimension smaller than n, i.\n");
        printf("e. the group is not an n-dimensional space group, the program\n");
        printf("will indicate an error, gives the rank, and exit.\n");
        printf("\n");
        printf("Options:\n");
        printf("-t    : The transforming matrix will be given as well.\n");
        printf("\n");
        printf("Cf. Extensions, Vector_systems\n");
        printf("Note: This program is a synonym for Extract -t.\n");
     }
     else{
     printf("Usage: %s 'file1' ['file2'] [-c] [-p] [-f] [-r [-D]] [-t=n]\n",argv[0]);
     printf("\n");
     printf("file1: bravais_TYP containing a space or (in case of option -r) a point group.\n");
     printf("file2: matrix_TYP only used with options -r and -t, cf. below.  \n");
     printf("\n");
     printf("Extracts the linear part of the affine (space) group in file1.\n");
     printf("Note this is the point group of the space group in case the space \n");
     printf("group is given in standard form.\n");
     printf("\n");
     printf("Options:\n");
     printf("-p    : extracts the linear part (default).\n");
     printf("-c    : extracts the translational part as a vector system (1-cozycle).\n");
     printf("-f    : do not calculate the formspace of the point group.\n");
     printf("-r    : reverses the process: Reads in the generators of the point group \n");
     printf("        of file1 and multiple vector systems for these generators from file2 \n");
     printf("        and outputs the resulting space groups. Overwrites any other option.\n");
     printf("-D    : this option only works together with '-r'. The spacegroups are written\n");
     printf("        to 'file.i' with i = 1, ...\n");
     printf("-t=n  : transform the space group with generators in  file1 into standard\n");
     printf("        form, ie. the translation lattice is transformed to Z^n. If a \n");
     printf("        parameter n>0 is specified, the transforming matrix will be\n");
     printf("        given as well.  NOTE: For the -t-option file2 must contain a \n");
     printf("        presentation for the point group.  If the translation lattice \n");
     printf("        is of smaller rank, it will give the rank. Synonymous command:\n");
     printf("        Standard_affine_form\n");
     printf("-T     : RESERVED\n");
     printf("\n");
     printf("Cf. Extensions, Vector_systems\n");
     }


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

  G = get_bravais(FILENAMES[0]);

  if (is_option('t') || type){

     if (FILEANZ == 2){
       XX = mget_mat(FILENAMES[1],&i);
       if (i>1){
          fprintf(stderr,"you should only give a single matrix as presention\n");
          exit(3);
       }
       X = XX[0];
       free(XX);

     }
     else{

       /* get the point group */
       H = init_bravais(G->dim-1);
       H->gen_no = G->gen_no;
       H->gen = (matrix_TYP **) malloc(H->gen_no * sizeof(matrix_TYP *));
       for (i=0;i<H->gen_no;i++){
          H->gen[i] = copy_mat(G->gen[i]);
          real_mat(H->gen[i],G->dim,G->dim-1);
          real_mat(H->gen[i],G->dim-1,G->dim-1);
       }


       /* calculate a presentation */
       base = get_base(H);

       strong = strong_generators(base,H,TRUE);

       X = pres(strong,H,NULL);

       for (i=0;i<H->dim;i++){
         free_mat(base[i]);
         free_bahn(strong[i]);
         free(strong[i]);
       }
       free(strong);
       free(base);
       free_bravais(H); H = NULL;
     }

     if (is_option('T')){
        real_mat(X,X->rows+G->dim-1,X->cols);
        for (i=0;i<G->dim-1;i++){
           X->array.SZ[X->rows-G->dim+1+i][0] = G->gen_no - i ;
        }
     }

     /* set the transformation matrix */
     T = translation_lattice(G->gen,G->gen_no,X);
     free_mat(X);
     X = NULL;

     /* echo the rank if needed */
     if (T->rows != T->cols){
        fprintf(stderr,"The rank of the translation lattice is %d\n",T->cols);
        exit(0);
     }

     real_mat(T,T->rows+1,T->cols+1);

     /* paranoia */
     rat2kgv(T);

     T->array.SZ[T->rows-1][T->cols-1] = T->kgv;
     Check_mat(T);

     if ((optionnumber('t')>0 && !type) ||
         (is_option('t') && type)){
        sprintf(comment,"transformation matrix for space group in %s",
                        FILENAMES[0]);
        put_mat(T,NULL,comment,2);
     }

     /* we should also calculate a positive definite, point group
        invariant form, and transform the group according to a mink_red
        of this form. left for further development */

     /* now transform the group */
     H = init_bravais(G->dim);
     H->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
     H->gen_no = G->gen_no;
     X = mat_inv(T);
     for (i=0;i<H->gen_no;i++){
        H->gen[i] = mat_kon(X,G->gen[i],T);
     }

     sprintf(comment,"space group of %s on Z^n",
                        FILENAMES[0]);
     put_bravais(H,NULL,NULL);
  }
  else{
     if (is_option('r')){
        COZ = mget_mat(FILENAMES[1], &cozanz);
        for (k = 0; k < cozanz; k++){
           rat2kgv(COZ[k]);
           Check_mat(COZ[k]);
           convert_cocycle_to_column(&COZ[k],1,G->dim,G->gen_no);

           /* is it a valid cozycle? */
           if ((G->dim * G->gen_no != COZ[k]->rows) || (COZ[k]->cols != 1)){
             fprintf(stderr,"The cozycle is not compatible to this point group\n");
             fprintf(stderr,"It should have %d * %d = %d rows\n",G->dim,G->gen_no,
                     G->dim*G->gen_no);
             exit(3);
           }
           H = init_bravais(G->dim+1);
           H->gen_no = G->gen_no;
           H->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));

           for (i=0;i<H->gen_no;i++){
              H->gen[i] = copy_mat(G->gen[i]);
              rat2kgv(H->gen[i]);
              Check_mat(H->gen[i]);
              real_mat(H->gen[i],H->dim,H->dim);
              iscal_mul(H->gen[i],COZ[k]->kgv);
              H->gen[i]->kgv = H->gen[i]->kgv * COZ[k]->kgv;
              for (j=0;j<H->dim-1;j++)
                 H->gen[i]->array.SZ[j][H->dim-1] = COZ[k]->array.SZ[i*(H->dim-1)+j][0];
              H->gen[i]->array.SZ[H->dim-1][H->dim-1] = COZ[k]->kgv;
              Check_mat(H->gen[i]);
           }
           if (is_option('D')){
              sprintf(file, "%s.%d", FILENAMES[0], k + 1);
              sprintf(comment, "space group to the point group of %s and the %d-th cozycle of %s",
                      FILENAMES[0], k+1, FILENAMES[1]);
              put_bravais(H, file, comment);
           }
           else{
              sprintf(comment, "space group to the point group of %s and the %d-th cozycle of %s",
                      FILENAMES[0], k+1, FILENAMES[1]);
              put_bravais(H, NULL, comment);
           }
           free_bravais(H);
           free_mat(COZ[k]);
        }
        H = NULL;
        free(COZ);
     }
     else if (is_option('p') || (!is_option('p') && ! is_option('c'))){
        /* extract the point group */
        H = init_bravais(G->dim-1);

        H->gen_no = G->gen_no;
        H->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
        for (i=0;i<H->gen_no;i++){
           H->gen[i] = copy_mat(G->gen[i]);
           real_mat(H->gen[i],H->dim,G->dim);
           real_mat(H->gen[i],H->dim,H->dim);
           rat2kgv(H->gen[i]);
           Check_mat(H->gen[i]);
        }

        if (!is_option('f'))
           H->form = formspace(H->gen,H->gen_no,1,&H->form_no);

        /* output */
        sprintf(comment,"point group to the space group of %s",FILENAMES[0]);
        put_bravais(H,NULL,comment);
        free_bravais(H);
        H = NULL;
     }
     if (is_option('c')){
        X = init_mat(G->gen_no * (G->dim-1),1,"");
        kgv = 1;
        for (i=0;i<G->gen_no;i++){
           rat2kgv(G->gen[i]);
           kgv = kgv*G->gen[i]->kgv/GGT(kgv,G->gen[i]->kgv);
        }


        /* set the cozycle */
        X->kgv = kgv;
        for(i=0;i<G->gen_no;i++)
           for (j=0;j<G->dim-1;j++)
              X->array.SZ[i*(G->dim-1)+j][0] =  X->kgv/G->gen[i]->kgv *
                                             G->gen[i]->array.SZ[j][G->dim-1];

        Check_mat(X);
        sprintf(comment,"cozycle to the group of %s",FILENAMES[0]);
        put_cocycle(X,G->dim-1,G->gen_no,NULL,comment);
     }
  }

  /* cleaning up */
  if (G!=NULL) free_bravais(G);
  if (H!=NULL) free_bravais(H);
  if (X!=NULL) free_mat(X);
  if (T!=NULL) free_mat(T);

  exit(0);

} /* main */
