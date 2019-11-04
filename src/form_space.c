#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int INFO_LEVEL;

int main (int argc, char *argv[])
{
  int i,anz, prime;
  int sym_opt;
  bravais_TYP *B;
  matrix_TYP **F;
  char comment[80];

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("Usage: Form_space 'file' [-a] [-s] [-p=prime]\n");
    printf("\n");
    printf("file: bravais_TYP containing the unimodular group G.\n");
    printf("\n");
    printf("Calculates a Z-basis for the space of matrices A with g^tr * A * g = A\n");
    printf("for all g in G. Default: A is symmetric.\n");
    printf("\n");
    printf("Options:\n");
    printf("-a:       all invariant matrices are calculated\n");
    printf("-s:       only the skew-symmetric invariant matrices are calculated\n");
    printf("-p=prime: the mod-p-invariant matrices, and a basis over Z/pZ is given.\n");
    printf("          The default prime is  101. \n");
    printf("\n");
    printf("Cf. Invar_space\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }
  B = get_bravais(FILENAMES[0]);
  sym_opt = 1;
  if(is_option('a') == TRUE)
    sym_opt = 0;
  if(is_option('s') == TRUE)
    sym_opt = -1;
  if(is_option('p'))
  {
     if(optionnumber('p') == 0)
        prime = 101;
     else
        prime = optionnumber('p');
     F = p_formspace(B->gen, B->gen_no, prime, sym_opt, &anz);
     sprintf(comment, "invariant martrix modulo %d", prime);
     printf("#%d\n", anz);
     for(i=0;i<anz;i++)
       put_mat(F[i], NULL, comment, 2);
  }
  else
  {
    F = formspace(B->gen, B->gen_no, sym_opt, &anz);
    printf("#%d\n", anz);
    for(i=0;i<anz;i++)
     put_mat(F[i], NULL, "invariant matrix", 2);
  }

  exit(0);
}
