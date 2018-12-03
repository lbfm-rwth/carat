#include <typedef.h>
#include <getput.h>
#include <gmp.h>
#include <zass.h>
#include <matrix.h>
#include <tools.h>
#include <bravais.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;

int main(int argc,char **argv){


  int i,
      j,
      k,
      anz,
      denominator,
    **words;

  bravais_TYP *P,
             **R;

  char comment[1000];

  matrix_TYP *coz,
            **RG;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) ||
     ((FILEANZ < 2))){
     printf("Usage: %s 'file1' 'file2' ['file3' .....]\n",argv[0]);
     printf("\n");
     printf("file1: bravais_TYP containing the point group G.\n");
     printf("file2: bravais_TYP containing the space group R_1.\n");
     printf("filex: bravais_TYP containing additional space groups.\n");
     printf("\n");
     printf("For each group in the files 'file2' 'file3' ... the program\n");
     printf("tries to find elements in the space group with the same\n");
     printf("linear part as the generators in file1.\n");
     printf("WARNING: If this is impossible, the program will run into an\n");
     printf("         infinite loop.\n");
     printf("\n");
     printf("CF. Extensions, Extract,Z_equiv.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  INFO_LEVEL = optionnumber('h');

  anz = FILEANZ-1;
  R = (bravais_TYP **) malloc( anz * sizeof(bravais_TYP *));

  /* read the files the user specified */
  P = get_bravais(FILENAMES[0]);
  for (i=0;i<anz;i++)
     R[i] = get_bravais(FILENAMES[i+1]);

  for (i=0;i<anz;i++){
     RG = (matrix_TYP **) malloc( R[i]->gen_no * sizeof(matrix_TYP*));
     words = (int **) malloc(P->gen_no * sizeof(int *));
     denominator = 1;
     for (j=0;j<R[i]->gen_no;j++){
       RG[j] = copy_mat(R[i]->gen[j]);
       RG[j]->cols--;
       RG[j]->rows--;
       Check_mat(RG[j]);
       if (!RG[j]->flags.Integral){
          fprintf(stderr,"The point group has to be integral\n");
          exit(3);
       }
       /* kgv2rat(R[i]->gen[j]); */
       rat2kgv(R[i]->gen[j]);
       denominator *= (R[i]->gen[j]->kgv / GGT(R[i]->gen[j]->kgv,denominator));
     }

     /* stick the rigth INTEGRAL cozycle at the end of the RG[j] */
     for (j=0;j<R[i]->gen_no;j++){
        RG[j]->cols++;
        RG[j]->rows++;
        for (k=0;k<RG[j]->rows-1;k++)
           RG[j]->array.SZ[k][R[i]->dim-1] = (denominator / R[i]->gen[j]->kgv) *
                   R[i]->gen[j]->array.SZ[k][R[i]->dim-1];
        RG[j]->array.SZ[R[i]->dim-1][R[i]->dim-1] = 1;
        Check_mat(RG[j]);
     }

     /* get the cozycle on the right generators */
     coz = reget_gen(RG,R[i]->gen_no,P,words,TRUE);

     /* the cozykle has to become the right denominator */
     coz->kgv = denominator;
     Check_mat(coz);

     /* output */
     sprintf(comment,"standartized cozycle for the group of %s",FILENAMES[i+1]);
     put_cocycle(coz,P->dim,P->gen_no,NULL,comment);


     for (j=0;j<R[i]->gen_no;j++){
        free_mat(RG[j]);
     }
     for (j=0;j<P->gen_no;j++){
        free(words[j]);
     }
     free(words);
     free(RG);
  } 

  /* clean up */
  free_bravais(P);
  for (i=0;i<anz;i++)
     free_bravais(R[i]);
  free(R);

 
  if (INFO_LEVEL & 12){
     SFLAG = 1;
  }

  if (INFO_LEVEL & 12){
     fprintf(stderr,"write pointer_statistics\n");
     pointer_statistics(0,0);
  }

  if (FALSE){
     put_bravais(P,NULL,NULL);
  }

  exit(0);

} /* main */
