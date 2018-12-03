#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <datei.h>
#include <longtools.h>
#include <bravais.h>
#include <idem.h>

int INFO_LEVEL;
extern int SFLAG;


int main(int argc,char **argv){

  bravais_TYP *G;

  matrix_TYP *form,
            **idem,
             *id;

  int i,
      anz,
      dim,
      dimc,
      dimcc,
      options[2];

  char comment[1000];

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s file [-d] [-o]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing generators of a finite unimodular group.\n");
     printf("\n");
     printf("Calculates the central primitive idempotents of the rational commuting (or \n");
     printf("equivalently envelloping) algebra of the group in file. (This is used\n");
     printf("for splitting the natural module Q^n into homogeneous components.)\n");
     printf("\n");
     printf("CAUTION: The program might give wrong answers for degrees bigger than 6.\n");
     printf("\n");
     printf("Options:\n");
     printf(" -d     : Outputs the dimensions of the commuting algebra and its centre.\n");
     printf(" -o     : Gives a Z-basis of the commuting algebra and its centre.\n");
     printf("\n");
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

  /* avoid silly mistakes */
  for (i=0;i<G->gen_no;i++){
     if (!G->gen[i]->flags.Integral){
        fprintf(stderr,"The generators for the bravais group have to be\n");
        fprintf(stderr,"integral.\n");
        exit(3);
     }
  }

  id = init_mat(G->dim,G->dim,"1");
  form = rform(G->gen,G->gen_no,id,101);

  idem = idempotente(G->gen,G->gen_no,form,&anz,&dimc,&dimcc,options);

  printf("#%d\n",anz);
  for (i=0;i<anz;i++){
     sprintf(comment,"%d-th idempotent for group in %s",i+1,FILENAMES[0]);
     put_mat(idem[i],NULL,comment,2);
  }

  if (is_option('o')){
     /* tell the user the centralizer */
     printf("#%d\n",dimc);
     for (i=0;i<dimc;i++){
        sprintf(comment,"centralizer of group in %s",FILENAMES[0]);
        put_mat(idem[i+anz],NULL,comment,2);
     }
     /* and it's centre */
     printf("#%d\n",dimcc);
     for (i=0;i<dimc;i++){
        sprintf(comment,"centre of centralizer for group in %s",FILENAMES[0]);
        put_mat(idem[i+anz+dimc],NULL,comment,2);
     }     
  }
  if (is_option('d')){
     printf("dimension of the centralizer: %d\n",dimc);
     printf("dimension of it's centre: %d\n",dimcc);
  }

  for (i=0;i<anz+dimc+dimcc;i++){
     free_mat(idem[i]);
  }

  if (idem != NULL) free(idem);
  free_bravais(G);
  free_mat(form);
  free_mat(id);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(0);

} /* main */
