#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <longtools.h>
#include <bravais.h>
#include <idem.h>
#include <voronoi.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;

int main(int argc,char **argv){

  bravais_TYP *G,
              *H;

  matrix_TYP *form,
             *id;

  int i,
      anz,
      almost,
      zclass,
      type,      /* TRUE iff the program is called via Bravais_type */
      options[2];

  char comment[1000],
       *symb;

  read_header(argc,argv);

  #define OTHER_NAME "Bravais_type"
  if (strlen(argv[0]) < strlen(OTHER_NAME)){
     /* it can't have been called via .....OTHER_NAME */
     type = FALSE;
  }
  else{
     type = (strcmp(argv[0]+(strlen(argv[0])-strlen(OTHER_NAME)),
             OTHER_NAME) == 0);
  }

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     if (!type){
        printf("Usage: %s 'file' [-i] [-B] [-t] [-f]\n",argv[0]);
     }
     else{
        printf("Usage: %s 'file' [-B] [-t] [-f]\n",argv[0]);
     }
     printf("\n");
     printf("file: bravais_TYP.\n");
     printf("\n");
     printf("Identify the family symbol of a finite unimodular group\n");
     printf("given in 'file'.\n");
     printf("\n");
     printf("Options:\n");
     if (!type){
     printf(" -i    : Identify the group in file even more,\n");
     printf("         ie. give the position and representative of it's\n");
     printf("         Bravais group in the catalog.\n");
     }
     printf(" -B    : Assume the group to be a Bravais group, ie. do not\n");
     printf("         calculate the bravais group of it. Note, the order of\n");
     printf("         the Bravais group must be correct.\n");
     if (!type){
     printf(" -t    : (Does only work in conjunction with -i) Give\n");
     printf("         the transformation matrix which conjugates\n");
     printf("         the given group to the group in the catalog.\n");
     }
     else{
     printf(" -t    : Give the transformation matrix which conjugates\n");
     printf("         the given group to the group in the catalog.\n");
     }
     printf(" -f    : recalculate the form space even if it is given.\n");
     if (!type){
     printf("         (has an effect only if given with the option -i)\n");
     }
     printf("\n");
     printf("Cf: Bravais_catalog, Datei\n");
     printf("Note: Bravais_type is a synonym for Symbol -i\n");
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

  if (G->form)
     long_rein_formspace(G->form,G->form_no,1);

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

  symb = symbol(G,form);

  printf("Symbol of the group: %s\n", symb);

  /* here it is important that G is really a bravais group */
  if (!is_option('B')){
     anz = 0;
     for (i=0;i<G->dim;i++)
        if (anz < form->array.SZ[i][i]) anz = form->array.SZ[i][i];

     if (anz > 10)
        H = bravais_group(G,TRUE);
     else
        H = bravais_group(G,FALSE);

     free_bravais(G);
     G = H;
     H = NULL;
  }

  if (is_option('i') || type){
     /* see if we need the formspace */
     if (is_option('f') || G->form_no == 0){
        for (i=0;i<G->form_no;i++) free_mat(G->form[i]);
        if (G->form_no>0) free(G->form);
        G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
     }

     free_mat(form);
     H = catalog_number(G,symb,&form,&almost,&zclass);
     printf("Position in the catalog:\n");
     printf("%d-th homogenously decomposable bravais group, %d-th Z-class\n",
             almost,zclass);
     put_bravais(H,NULL,NULL);
     free_bravais(H);

     if (is_option('t')) put_mat(form,NULL,"transformation matrix",2);
  }

  free(symb);
  free_bravais(G);
  free_mat(form);
  free_mat(id);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(0);

} /* main */
