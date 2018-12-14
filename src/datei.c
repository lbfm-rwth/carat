#include "typedef.h"
#include "tools.h"
#include "datei.h"
#include "getput.h"
#include "matrix.h"

int main (int argc, char *argv[])

{
  char *filename;
  char string[80];
  char *fn,
       *family_name;
  symbol_out **B;
  bravais_TYP *S;
  matrix_TYP *X;
  int ad_no = 0;   /* number of homogenously decomposable groups in this family */
  FILE *outfile;
  int i,j,c;
  
  extern bravais_TYP *Z_class();
  extern symbol_out *read_symbol();
  extern symbol_out *get_symbol();

  /*  scan_argv (argc, argv, &filename);  */
  /* changed to conform with the rest of the carat package */
  read_header(argc,argv);
  if (FILEANZ > 0)
     filename = FILENAMES[0];
  else
     filename = NULL;

  if ((is_option('h') && optionnumber('h') == 0) ||
      FILEANZ > 1){
     printf("Usage: %s [file]\n",argv[0]);
     printf("\n");
     printf("file: (OPTIONAL) contains a set of commands for Datei, which are otherwise\n");
     printf("      asked from stdin.\n");
     printf("\n");
     printf("Accesses the catalog of bravais groups up to dimension 6.\n");
     printf("\n");
     printf("The first input for the catalog is the family symbol, whose grammar is\n");
     printf("described in detail in the CARAT manual. It is build up from atomic symbols\n");
     printf("which are seperated by `,' or `;' to indicate diagonals and direct products.\n");
     printf("We just state the two important rules here:\n");
     printf("\n");
     printf(" LIST OF 'ATOMS':\n");
     printf(" dim1: 1\n");
     printf(" dim2: 2-1  2-1'  2-2  2-2'\n");
     printf(" dim3: 3\n");
     printf(" dim4: 4-1  4-1'  4-2  4-2'  4-3  4-3'\n");
     printf(" dim5: 5-1  5-2\n");
     printf(" dim6: 6-1  6-2  6-2'  6-3  6-3'  6-4  6-4'\n");
     printf("\n");
     printf("Meaning of `,':                        (X 0 0)\n");
     printf(" A,A,A stands for groups of the form   (0 X 0) with X in A\n");
     printf("                                       (0 0 X)\n");
     printf("Meaning of `;':\n");
     printf(" A;B stands for groups of the form (X 0)\n");
     printf("                                   (0 Y) with X in A and Y in B.\n");
     printf("NOTE: The primed atoms 2-1' and 2-2' only occur in multiples, seperated by `,'.\n");
     printf("Examples: 2-1,2-1;2-1     (degree 6)\n");
     printf("          1;1;1           (degree 3)\n");
     printf("          3,3             (degree 6)\n");
     printf("          4-1;1           (degree 5)\n");
     printf("\n");
     printf("Cf.: Symbol, Bravais_type.\n");
     printf("Note: Bravais_catalog is a synonym for Datei.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }


  B = (symbol_out **) malloc(1 *sizeof(symbol_out *));
  B[0] = read_symbol(filename);
  ad_no++;
  family_name = (char *) malloc(80 *sizeof(char));
  strcpy(family_name, B[0]->fn);
  /*
  fprintf(stderr,"B[0]->fn: %s\n",B[0]->fn);
  */
  family_name = strstr(family_name, "dim");
  family_name = family_name+5;
  get_zentr(B[0]);
  while(B[ad_no-1]->fn != NULL)
  {
    B = (symbol_out **) realloc(B, (ad_no+1) *sizeof(symbol_out *));
    filename = B[ad_no-1]->fn;
    B[ad_no] = get_symbol(filename);
    ad_no++;
  }
  printf("The crystal-family %s contains %d homogeneously decomposable bravais-groups with\n", family_name, ad_no);
  j=0;
  for(i=0; i<ad_no; i++)
  {
    if(j != 0)
    {
      if(j == ad_no -1)
         printf(" resp. ");
      else
        printf(", ");
    }
      printf(" %d", (B[i]->grp->zentr_no+1));
      j = 1;
  }
  printf("\n Z-classes of Bravais groups.\n");

  
  printf("Do you want to calculate bravais-groups? (y or n): ");
  fn = (char *) malloc(80 *sizeof(char));
  c=scanf( "%[ \t\n]", fn);
  c=scanf( "%[^\n]", fn);
  while(strncmp(fn, "y", 1) != 0 && strncmp(fn, "n", 1) != 0)
  {
    c=scanf( "%[ \t\n]", fn);
    c=scanf( "%[^\n]", fn);
  }
  if(strncmp(fn, "n", 1) == 0)
    exit(0);

    /*-----------------------------------------------------*\
    | read and open output-file                             |
    \*-----------------------------------------------------*/
  printf("Please input a filename (stdout = standard output): ");
  c=scanf( "%[ \t\n]", fn);
  c=scanf( "%[^\n]", fn);
  if(strncmp(fn, "stdout", 6) == 0)
  fn = NULL;
  if(fn == NULL)
    outfile = stdout;
  else
    outfile = fopen( fn, "w");

  printf("Which bravais-groups should be printed? (a(ll) or s(election): ");
  c=scanf( "%[ \t\n]", string);
  c=scanf( "%[^\n]", string);
  while(strncmp(string, "a", 1) != 0 && strncmp(string, "s", 1) != 0)
  {
    c=scanf( "%[ \t\n]", string);
    c=scanf( "%[^\n]", string);
  }

  if(strncmp(string, "a", 1) == 0)
  {
  for(i=0; i<ad_no; i++)
  {
     fput_bravais(outfile, B[i]->grp, "homogenously decomposable bravais-group");
     fflush(outfile);
     for(j=0; j<B[i]->grp->zentr_no; j++)
     {
        X = B[i]->grp->zentr[j];
        S = Z_class(B[i]->grp, X);
        fput_bravais(outfile, S, "not homogenously decomposable bravais-group");
        fflush(outfile);
        free_bravais(S);
     }
  }
  exit(0);
  }

  strncpy(string, "y", 1);
  while(strncmp(string, "y", 1) == 0)
  {
     printf("Please enter index i of homogeneously decomposable bravais-group, 1<= i<= %d: ", ad_no);
     c=scanf("%d", &i);
     i--;
     if(i<0 || i>= ad_no)
        printf("There is no homogenously decomposable bravais-group of this index\n");
     else
     {
       printf("Please enter index j of bravais-group belonging to this homogenously decomposable\n");
       printf("1 <= j =< %d: ", (B[i]->grp->zentr_no+1));
       c=scanf("%d", &j);
       j--;
       if(j== 0)
        fput_bravais(outfile, B[i]->grp, "homogeneously decomposable bravais-group");
       else
       {
         j--;
         if(j<0 || j>=B[i]->grp->zentr_no)
             printf("There is no bravais-group of this index\n");
         else
         {
            X = B[i]->grp->zentr[j];
            S = Z_class(B[i]->grp, X);
            fput_bravais(outfile, S, "not homogeneously decomposable bravais-group");
            free_bravais(S);
         }
       }
     }
     printf("Do you want further bravais-groups? (y or n): ");
     c=scanf( "%[ \t\n]", string);
     c=scanf( "%[^\n]", string);
     while(strncmp(string, "y", 1) != 0 && strncmp(string, "n", 1) != 0)
     {
       c=scanf( "%[ \t\n]", string);
       c=scanf( "%[^\n]", string);
     }

     /* inserted to enable multiple output files, tilman 05/06/97 */
     if (strncmp(string,"y",1) == 0){
        printf("If you want a different file, insert it: (no/filename)");
	if (fn == NULL) fn = (char *) malloc(1024*sizeof(char));
        c=scanf( "%[ \t\n]", fn);
        c=scanf( "%[^\n]", fn);
        if (strcmp(fn,"no") != 0){
	   fclose (outfile);
           outfile = fopen( fn, "w");
        }
     }
  }

if ( outfile != stdout ) {
	fclose (outfile);
	}

exit(0);
}
/*{{{}}}*/
