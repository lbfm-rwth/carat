#include <ZZ.h>
#include <typedef.h>
#include <bravais.h>
#include <base.h>

int INFO_LEVEL;
extern int SFLAG;
extern int IDEM_NO;
boolean GRAPH = FALSE;
boolean GRAPH_DEBUG = FALSE;

int main(int argc,char **argv){

  bravais_TYP *G,
             **Classes;

  bahn **ST;

  int i,
      no,
      adnumber = 0,
      second_number = 0;

  char comment[1000],
       file[1000];

  matrix_TYP **basis;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s 'file' [-D] [-a] [-f]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a group. The order of the group\n");
     printf("      must be specified.\n");
     printf("\n");
     printf("Splits the Q-class of the group in 'file' into Z-classes and\n");
     printf("gives a representative for each Z-class.\n");
     printf("\n");
     printf("Options:\n");
     printf("-a    : Only calculates representatives of the homogeneously \n");
     printf("        decomposable Z-classes and basis transformations \n");
     printf("        to the remaining Z-classes (as centerings).\n");
     printf("        This option is faster.\n");
     printf("-D    : (For debugging:) Output is  split into  different\n");
     printf("        files 'file'.1, 'file'.2, ... matching the\n");
     printf("        representative groups.\n");
     printf("-f    : recalculate the formspace even if it is given.\n");
     printf("-q    : quiet mode.\n");
     printf("\n");
     printf("WARNING: This program might be very time consuming, especially\n");
     printf("         if the group has many Z-classes. \n");
     printf("\n");
     printf("Cf. Order, Is_finite.\n");
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

  if (G->form == NULL ||
      G->form_no == 0 ||
      is_option('f')){
     G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
  }

  if (G->order == 0){
    basis = get_base(G);
    ST = strong_generators(basis,G,FALSE);
    G->order = 1;
    for (i=0;i<G->dim;i++){
      G->order *= ST[i]->length;
      free_mat(basis[i]);
      free_bahn(ST[i]);
      free(ST[i]);
    }
    factorize_new(G->order,G->divisors);
    free(ST);
    free(basis);
  }

  Classes = q2z(G,&no,is_option('a'), NULL, is_option('q'));

  fprintf(stderr, "####### There are %d classes of groups\n", no);
  for (i=0;i<no;i++){
     if (is_option('a') || Classes[i + no]){
        adnumber++;
        second_number = 0;
     }
     second_number++;
         
     sprintf(comment,"%d-th homogenously dec. group, %d zclass",
                       adnumber,second_number);
     if (is_option('D')){
        sprintf(file,"%s.%d.%d",FILENAMES[0],adnumber,second_number);
        put_bravais(Classes[i],file,comment);
     }
     else{
        put_bravais(Classes[i],NULL,comment);
     }
  }

  free_bravais(G);
  for (i=0;i<no;i++){
     free_bravais(Classes[i]);
  }
  free(Classes);
  cleanup_prime();

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(0);

} /* main */
