#include <ZZ.h>
#include <typedef.h>
#include <bravais.h>

int INFO_LEVEL;
extern int SFLAG;
extern int IDEM_NO;

main(int argc,char **argv){

  bravais_TYP *G,
             **Classes;

  int i,
      no;

  char comment[1000],
       file[1000];


  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s 'file' [-D] [-a] [-f]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a group. Just must specify it's order.\n");
     printf("\n");
     printf("Calculates representatives of the Z-classes of groups in\n");
     printf("the same Q-class as the group in file.\n");
     printf("\n");
     printf("Options:\n");
     printf("-a    : just calculate the homogeneously decomposable groups and\n");
     printf("        centerings, which which rise to the Z-classes.\n");
     printf("        This option is faster.\n");
     printf("-D    : (meant for debugging) write the output on different\n");
     printf("        files GROUP.1, GROUP.2, ... for each group.\n");
     printf("-f    : recalculate the formspace even if it is given.\n");
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

  Classes = q2z(G,&no,is_option('a'));

  printf("####### There are %d classes of groups\n",no);
  for (i=0;i<no;i++){
     sprintf(comment,"%d-th zclass",i+1);
     if (is_option('D')){
        sprintf(file,"GROUP.%d",i+1);
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
