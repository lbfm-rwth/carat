#include <typedef.h>
#include <getput.h>
#include <bravais.h>
#include <matrix.h>
#include <base.h>
#include <datei.h>
#include <idem.h>
#include <longtools.h>

int INFO_LEVEL;
extern int SFLAG;

void main(int argc,char **argv){

  bravais_TYP *G = NULL,
              *H;

  matrix_TYP *F,
             *FI,
             *ID,
             *T= NULL;

  lattice_element **RES;

  int almost,
      zclass,
      number,
      i,
      j,
      no_N_classes;

  char comment[1000],
      *symb = NULL;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) ||
     ((FILEANZ != 1))){
     printf("Usage: %s 'file' [-S]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a Bravais group G.\n");
     printf("\n");
     printf("Lists representatives of the Bravais subgroups (resp. Bravais supergroups \n");
     printf("with option -S) of the BRAVAIS GROUP G under the conjugation action of\n");
     printf("the normalizer of G in GL_n(Z).\n");
     printf("Remark: Without option -G only the symbols of the Bravais groups are given.\n");
     printf("\n");
     printf("Options:\n");
     printf("-S:   Search for all supergroups\n");
     printf("-G:   Calculate generators for the groups as well\n");
     printf("\n");
     printf("Cf. Bravais_catalog, Symbol.\n");
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

  if (FILEANZ == 1){
     G = get_bravais(FILENAMES[0]);

     /* cause the program to recalculate the formspace if desired */
     if (is_option('f')){
        for (i=0;i<G->form_no;i++){
           free_mat(G->form[i]);
        }
        free(G->form);
        G->form=NULL;
        G->form_no = 0;
     }

     /* we desperately will need the space of forms */
     if (G->form_no == 0 || G->form == NULL){
        G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
     }

     /* if (!is_option('b')){
        H = bravais_group(G,FALSE);
        free_bravais(G);
        G = H;
        H = NULL;
     } */

     ID = init_mat(G->dim,G->dim,"1");
     F = rform(G->gen,G->gen_no,ID,101);
     symb = symbol(G,F);
     H = catalog_number(G,symb,&T,&almost,&zclass);
     free_bravais(H);
     free_mat(F);
     free_mat(ID);
  }

  if (is_option('S')){
     RES = super_lattice(symb,G->dim,almost,zclass,&number,is_option('G'));

     if (!is_option('G')){
        printf("Bravais groups which contain a Z-equivalent subgroup\n");
        for (i=0;i<number;i++){
           printf("Symbol: %s  homogeneously d.: %d zclass: %d\n",
                  RES[i]->symbol,RES[i]->almost,RES[i]->zclass);
        }
     }
     else{
        /* we now have to care for the transformation matrix */
        printf("bravais supergroups of the groups G in %s\n",FILENAMES[0]);
        printf("number of Z-conjugacy classes: %d\n",number);
        no_N_classes = 0;
        for (i=0;i<number;i++){
           no_N_classes += RES[i]->N_orbits;
        }
        printf("number of N_GLn(Z)(G)-conjugacy classes: %d\n",no_N_classes);
        for (i=0;i<number;i++){
           for (j=0;j<RES[i]->N_orbits;j++){
              F = mat_mul(T,RES[i]->TR[j]);
              FI = long_mat_inv(F);
              H = konj_bravais(RES[i]->grp,FI);
              free_mat(F);
              free_mat(FI);
              sprintf(comment,"%s %d %d, %d-th N-orbit",RES[i]->symbol,
                           RES[i]->almost,RES[i]->zclass,j+1);


              if (is_option('A')){
                 sprintf(comment,"%s.%s.%d.%d.%d",FILENAMES[0],RES[i]->symbol,
                           RES[i]->almost,RES[i]->zclass,j+1);
                 put_bravais(H,comment,NULL);
              }
              else{
              put_bravais(H,NULL,comment);
              }
              free_bravais(H);
           }
        }
     }

  }
  else{
     RES = lattice(symb,G->dim,almost,zclass,&number,is_option('G'));

     if (!is_option('G')){
        printf("Bravais subgroups of the groups G in %s\n",FILENAMES[0]);
        for (i=0;i<number;i++){
           printf("Symbol: %s  homogeneously d.: %d zclass: %d\n",
                  RES[i]->symbol,RES[i]->almost,RES[i]->zclass);
           printf("Number of occurences in the group: %d\n",RES[i]->alpha);
           printf("Number of orbits under the normalizer: %d\n",
                                              RES[i]->N_orbits);
        }
     }
     else{
        /* we now have to care for the transformation matrix */
        printf("bravais subgroups of the groups G in %s\n",FILENAMES[0]);
        printf("number of Z-conjugacy classes: %d\n",number);
        no_N_classes = 0;
        for (i=0;i<number;i++){
           no_N_classes += RES[i]->N_orbits;
        }
        printf("number of N_GLn(Z)(G)-conjugacy classes: %d\n",no_N_classes);
        F = mat_inv(T);
        free_mat(T);
        T = F;
        for (i=0;i<number;i++){
           for (j=0;j<RES[i]->N_orbits;j++){
              F = mat_inv(RES[i]->TR[j]);
              ID = mat_mul(T,F);
              free_mat(F);
              H = konj_bravais(RES[i]->grp,ID);
              free_mat(ID);
              sprintf(comment,"%s %d %d, %d-th N-orbit",RES[i]->symbol,
                           RES[i]->almost,RES[i]->zclass,j+1);
              put_bravais(H,NULL,comment);
              free_bravais(H);
           }
        }
     }
  }


  /* clean up the lattice */
  for (i=0;i<number;i++){
     free_lattice_element(RES[i]);
  }
  free(RES);

  /* clean up */
  if (G != NULL) free_bravais(G);
  if (T != NULL) free_mat(T);
  if (symb != NULL) free(symb);

  if (INFO_LEVEL & 12){
     fprintf(stderr,"write pointer_statistics\n");
     pointer_statistics(0,0);
  }

  exit(0);

} /* main */
