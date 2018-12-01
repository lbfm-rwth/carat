#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int i,
      j,
      anz,
      epsilon,
      expected_dimension;

  int sym_opt;
  bravais_TYP *B;
  matrix_TYP **F;
  char comment[80];

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
     printf("Usage: %s 'file' [-a] [-s] [-e=n] [-d=n]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a finite unimodular group G.\n");
     printf("\n");
     printf("Calculates a Z-basis for the space of matrices with g^{tr}Ag = A for all\n");
     printf("g in G by a quick seminumerical algorithm with some random features. By \n");
     printf("default all symmetric invariant matrices are calculated,\n");
     printf("\n");
     printf("Options:\n");
     printf("-a:      All invariant matrices are calculated.\n");
     printf("-s:      All skew-symmetric invariant matrices are calculated.\n");
     printf("-e=n:    Changes the default value of control parameter from 100 to n>100. \n");
     printf("-d=n:    The program will calculate n+1 elements F_0,..., F_n of the \n");
     printf("         invariant space, then it calculates a Z-basis of the Q-span of\n");
     printf("         F_0, ..., F_n intersected with Z^nxn.\n");
     printf("         This option is only usefull if one knows the dimension of the \n");
     printf("         space of invariant forms. By default this number is calculated \n");
     printf("         by modular arithmetic.\n");
     printf("\n");
     printf("Cf. Form_space.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  B = get_bravais(FILENAMES[0]);
  sym_opt = 1;

  if(is_option('a') == TRUE){
    sym_opt = 0;
  }
  if(is_option('s') == TRUE){
    sym_opt = -1;
  }

  if (is_option('d')){
     expected_dimension = optionnumber('d');
  }
  else{
     /* calculating the dimension over a prime field */
     F = p_formspace(B->gen,B->gen_no,1949,sym_opt,&expected_dimension);
  }

  if(optionnumber('e') == 0)
     epsilon = 101;
  else
     epsilon = optionnumber('e');

  F = invar_space(B->gen, B->gen_no,expected_dimension,sym_opt,epsilon,&anz);
  printf("#%d\n", anz);
  for(i=0;i<anz;i++)
  put_mat(F[i], NULL, "invariant matrix", 2);

  exit(0);
}
