#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"
#include"datei.h"


/* ---------------------------------------- */
/* Convert a single matrix in MAPLE format. */
/* ---------------------------------------- */
static void conv_mat_to_maple(matrix_TYP *A, 
                              int i,
                              const char *name){
    int j, k;


    Check_mat(A);
    if (!A->flags.Integral){
       rat2kgv(A);
    }
    printf("%s%d := array(1..%d,1..%d,[\n", name, i+1, A->rows, A->cols);
    for (j=0;j<A->rows;j++){
       printf("[");
       for (k=0;k<A->cols;k++){
          printf(" %d",A->array.SZ[j][k]);
          if (!A->flags.Integral){
             printf("/%d",A->kgv);
          }
          if (k!=(A->cols-1)){
             printf(",");
          }
       }
       if (j != (A->rows-1)){
          printf("],\n");
       }
       else{
          printf("]");
       }
    }
    printf("]):\n");
}



/* ---------------------------------------- */
/* Convert a single matrix in GAP format.   */
/* ---------------------------------------- */
static void conv_mat_to_gap(matrix_TYP *A,
                            int i,
                            const char *name){
   int j, k;


    Check_mat(A);
    if (!A->flags.Integral){
       rat2kgv(A);
    }
    printf("%s%d := [", name, i+1);
    for (j=0;j<A->rows;j++){
       printf("[");
       for (k=0;k<A->cols;k++){
          printf(" %d",A->array.SZ[j][k]);
          if (!A->flags.Integral){
             printf("/%d",A->kgv);
          }
          if (k!=(A->cols-1)){
             printf(",");
          }
       }
       if (j != (A->rows-1)){
          printf("],\n");
       }
       else{
          printf("]");
       }
    }
    printf("];;\n");
}



/* ---------------------------------------- */
/* Convert a single matrix in TeX format.   */
/* ---------------------------------------- */
static void conv_mat_to_tex(matrix_TYP *A){

   int j, k;


    Check_mat(A);
    if (!A->flags.Integral){
       rat2kgv(A);
    }
    printf("\\left(\\begin{array}{");
    for (j = 0; j < A->cols; j++)
       printf("c");
    printf("}\n");
    for (j = 0; j < A->rows; j++){
       for (k = 0;k < A->cols; k++){
          if (A->array.SZ[j][k] % A->kgv == 0)
             printf("%d ", A->array.SZ[j][k]/ A->kgv);
	  else
	     printf("%d/%d ", A->array.SZ[j][k], A->kgv);
	  if (k != A->cols - 1)
	     printf("& ");
       }
       printf("\\\\\n");
    }
    printf("\\end{array}\\right)\n");
}



/* ---------------------------------------- */
/* Main function.                           */
/* ---------------------------------------- */
int main (int argc, char *argv[])
{
  int i,
      anz,
      flag;

  matrix_TYP **A;

  bravais_TYP *G;

  const char *name;


  read_header(argc, argv);
  if(FILEANZ < 1 || (is_option('h') && optionnumber('h') == 0) )
  {
    printf("Usage: Conv 'file' ['name'] -g -G -m -M -t -T\n");
    printf("\n");
    printf("file: matrix_TYP or bravais_TYP (only generators are converted).\n");
    printf("name: (optional) a word\n");
    printf("\n");
    printf("Converts the matrix_TYP into the format of other\n");
    printf("programs available. The options specify the program,\n");
    printf("and only one of them is allowed at one time.\n");
    printf("First of all the number of matrices is printed.\n");
    printf("The matrices are called 'name.1', 'name.2', ... or\n");
    printf("if name isn't given 'a.1', 'a.2', ...\n");
    printf("\n");
    printf("-g or -G: GAP-format\n");
    printf("-m or -M: MAPLE-format\n");
    printf("-t or -T: TeX-format (the matrices aren't named and\n");
    printf("          the number of matrices isn't given)\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  G = get_bravais(FILENAMES[0]);
  anz = G->gen_no;
  A = G->gen;

  if (FILEANZ < 2)
     name = "a";
  else
     name = FILENAMES[1];

  flag = 0;
  if (is_option('t') || is_option('T'))
     flag = 3;
  if (is_option('m') || is_option('M'))
     flag = 2;
  if (is_option('g') || is_option('G'))
     flag = 1;

  if (flag > 0 && flag < 3)
     printf("%sNo := %i", name, anz);

  switch(flag){
     case 2:
        printf(":\n");
        for (i = 0; i < anz; i++)
           conv_mat_to_maple(A[i], i, name);
        break;

     case 1:
        printf(";;\n");
        for (i = 0; i < anz; i++)
           conv_mat_to_gap(A[i], i, name);
        break;

     case 3:
        for (i = 0; i < anz; i++)
	   conv_mat_to_tex(A[i]);
        break;

     default:
        printf("\nYou must specify at least one option. Try -h first.\n");
  }

  A = NULL;
  free_bravais(G);

  exit(0);
}










