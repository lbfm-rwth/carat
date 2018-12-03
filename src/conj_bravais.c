#include"typedef.h"
#include"getput.h"
#include"matrix.h"
#include"bravais.h"

int INFO_LEVEL;

int main (int argc, char *argv[])
{

   matrix_TYP **T,
               *X;
   bravais_TYP *G, *B;
   int anz;


   read_header(argc, argv);
   if(FILEANZ != 2)
   {
      printf("Usage: Conj_bravais 'file1' 'file2' [-i] \n");
      printf("\n");
      printf("file1: bravais_TYP containing the group G.\n");
      printf("file2: matrix_TYP containing a single matrix T.\n");
      printf("\n");
      printf("The program calculates TGT^{-1} and transforms all relevant\n");
      printf("data of G. The result is a bravais_TYP, that is written to\n");
      printf("standart output.\n");
      printf("\n");
      printf("Options:\n");
      printf(" -i     : Transform the group by T^{-1} G T.\n");
      printf("\n");
      if (is_option('h')){
         exit(0);
      }
      else{
         exit(31);
      }
   }
   G = get_bravais(FILENAMES[0]);
   T = mget_mat (FILENAMES[1], &anz);
   rat2kgv(T[0]);

   if (is_option('i')){
      X = mat_inv(T[0]);
      free_mat(T[0]);
      T[0] = X;
   }

   B = konj_bravais(G, T[0]);
   put_bravais(B, NULL, "konjugated bravais group");

   exit(0);
}
