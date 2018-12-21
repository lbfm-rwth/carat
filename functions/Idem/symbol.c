#include"typedef.h"
#include"idem.h"
#include"longtools.h"
#include"getput.h"
#include"voronoi.h"
#include"bravais.h"
#include"datei.h"
#include"matrix.h"
#include"orbit.h"
#include"tools.h"

extern int INFO_LEVEL;
/*****************************************************************************
@
@-----------------------------------------------------------------------------
@ FILE: symbol.c
@-----------------------------------------------------------------------------
@
*****************************************************************************/

/*****************************************************************************
@
@-----------------------------------------------------------------------------
@
@ static constituent *homogenous(bravais_TYP *G,matrix_TYP *F,int *anz)
@
@-----------------------------------------------------------------------------
@
*****************************************************************************/
static constituent *homogenous(bravais_TYP *G,matrix_TYP *F,int *anz)
{
   int i,
       j,
       rang,
       dimc,
       dimcc;

   matrix_TYP **idem,
              **solution,
               *colspace,
               *tmp;

   constituent *result;

   /* calculate the central primitive idempotents */
   idem = idempotente(G->gen,G->gen_no,F,anz,&dimc,&dimcc,FALSE);

   /* now we know allmost everything, just calculate a bit */
   result = (constituent *)malloc(anz[0] * sizeof(constituent));

   if (anz[0] > 1){
      for (i=0;i<anz[0];i++){

         /* calculate the dimension of the module the constiuent is acting on,
            and a basis of it */
         colspace = tr_pose(idem[i]);
         rang = tgauss(colspace);
         result[i].group = init_bravais(rang);
         real_mat(colspace,rang,colspace->cols);
         tmp = tr_pose(colspace);
         free_mat(colspace);
         colspace = tmp;
         tmp = NULL;
          
         result[i].group->gen = (matrix_TYP **) malloc(G->gen_no *
                                  sizeof(matrix_TYP *));
         result[i].group->gen_no = G->gen_no;

         /* transform the generators */
         for (j=0;j<G->gen_no;j++){
            tmp = mat_mul(G->gen[j],colspace);
            solution = long_solve_mat(tmp,colspace);

            if (solution[1] != NULL){
               fprintf(stderr,"error in homogenous: there shoudn't be\n");
               fprintf(stderr,"an homogenous solutio\n");
               exit(3);
            }
            result[i].group->gen[j] = solution[0];
            free(solution);
            free_mat(tmp);
         }

         /* now deal with the rest (not implemented yet) */
         result[i].centralizer = NULL;
         result[i].dimc=0;
         result[i].ccentralizer = NULL;
         result[i].dimcc=0;

         /* and calculate the formspace */
         result[i].group->form = formspace(result[i].group->gen,
                          result[i].group->gen_no,1,&result[i].group->form_no);

         free_mat(colspace);
      }

      /* clean up the result returned by idem */
      for (i=0;i<anz[0]+dimc+dimcc;i++) free_mat(idem[i]);
   }
   else{
      /* just copy the result */
      result[0].group = copy_bravais(G);

      result[0].dimc = dimc;
      result[0].dimcc = dimcc;
      result[0].centralizer = (matrix_TYP **) malloc(dimc
                                              * sizeof(matrix_TYP*));
      result[0].ccentralizer = (matrix_TYP **) malloc(dimcc
                                              * sizeof(matrix_TYP*));
      for (i=0;i<dimc;i++)
        result[0].centralizer[i] = idem[i+1];
      for (i=0;i<dimcc;i++)
        result[0].ccentralizer[i] = idem[i+1+dimc];

      /* we promised the user that there will be a formspace */
      if (G->form_no ==0)
         result[0].group->form = formspace(result[0].group->gen,
                          result[0].group->gen_no,1,&result[0].group->form_no);

      free_mat(idem[0]);
   }

   free(idem);

   return result;

}

/*****************************************************************************
@
@-----------------------------------------------------------------------------
@
@ static char *identify_hom(bravais_TYP *G,int clear)
@
@-----------------------------------------------------------------------------
@
*****************************************************************************/
static char *identify_hom(bravais_TYP *G,int clear)
{
  int i, c,
       found=0,
       flag,
       multiplicity,
      *list;

   static int atom_no;

   FILE *atom_file;

   char *result,
         tmp[20],
        *pp,
         bravais_file[1024], format[1024];

   static symbol_out *Atoms;

   /* this program does only work for groups up to dimension 6
      (it may even give wrong results otherwise), so tell the user */
   if (G->dim > 6){
      fprintf(stderr,"Error (in identify_hom), this program does only work\n");
      fprintf(stderr,"in dimension up to 6.\n");
      exit(3);
   }

   if (atom_no == 0){
       /* read the file with all atoms */
       get_data_dir(bravais_file, "tables/symbol/atom_list");
       atom_file = fopen(bravais_file,"r");

       if (atom_file == NULL){
	 fprintf(stderr, "I didn't find my library file %s. Exiting.\n",
		 bravais_file);
          exit(4);
       }

       c=fscanf(atom_file,"%d\n",&atom_no);

       /* now read the symbols */
       Atoms = (symbol_out *) malloc(atom_no * sizeof(symbol_out));
       for (i=0;i<atom_no;i++){
          Atoms[i].fn = (char *) calloc(20,sizeof(char));
          c=fscanf(atom_file,"%s\n",Atoms[i].fn);
       }
       fclose(atom_file);

       /* and their respective groups */
       get_data_dir(format, "tables/symbol/%s");
       for (i=0;i<atom_no;i++){
          sprintf(bravais_file, format, Atoms[i].fn);
          Atoms[i].grp = get_bravais(bravais_file);
          long_rein_formspace(G->form,G->form_no,1);
       }

       if (INFO_LEVEL & 4){
          for (i=0;i<atom_no;i++){
             printf("%s\n",Atoms[i].fn);
             put_bravais(Atoms[i].grp,NULL,NULL);
          }
       }
   }

   if (G->order == 0){
       fprintf(stderr,"this feature hasn't been implemented yet,\n");
       fprintf(stderr,"must specify a group order in identify_hom\n");
       exit(3);
   }

   list = (int *) malloc(atom_no * sizeof(int));

   /* check up the cheap stuff, ie. order, dimension ... */
   /* bare in mind that this is sufficient up to dimension 6 */
   for (i=0;i<atom_no;i++){
      flag = (Atoms[i].grp->order == G->order);
      flag = flag && ((G->dim % Atoms[i].grp->dim)==0);
      if (flag){
         list[found] = i;
         found++;
      }
   }

   /* just look whether these bravais groups really do give
      different families, or they look like (ie look for a,b,c's..)*/
   if (found > 1){
      /* to be implemented, it doesn't matter for dimensions up to 6 */
   }

   /* now list contains a list of representatives which are possible.
      now rule out the rest */
   if (found == 1){
      multiplicity = G->dim / Atoms[list[0]].grp->dim;
   }
   else{
      fprintf(stderr,"this feature hasn't been implemeneted yet\n");
      fprintf(stdout,"please report this to carat@momo.math.rwth-aachen.de\n");
      fprintf(stdout,"together with the output:\n");
      put_bravais(G,NULL,"ERROR IN: identify_hom");
      exit(3);
   }

   /* now we know the Q-Class of this homogenous module */
   result = (char *) malloc(20*multiplicity*sizeof(char));
   /* bare in mind to remove a,b,c and stuff from the end of the type */
   pp = strchr(Atoms[list[0]].fn,'a');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'b');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'c');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'d');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'e');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'f');
   if (pp == NULL) pp = strchr(Atoms[list[0]].fn,'g');
   if (pp == NULL){
       sprintf(result,"%s",Atoms[list[0]].fn);
   }
   else{
       strncpy(result,Atoms[list[0]].fn,pp - Atoms[list[0]].fn);
       result[pp-Atoms[list[0]].fn] = 0;
   }
   sprintf(tmp,"%s",result);
   for (i=1;i<multiplicity;i++){
      sprintf(result,"%s,%s",result,tmp);
   }

   /* if we got the order, clean up the static allocated memory */
   if (clear){
      for (i=0;i<atom_no;i++){
         free_bravais(Atoms[i].grp);
         free(Atoms[i].fn);
      }
      free(Atoms);
      atom_no = 0;
      Atoms = NULL;
   }

   free(list);

   return result;
}

static void bubblesort(char **S,int n)
{
  int i=1;

  char *tmp;

  while (i<n){
     if (strcmp(S[i-1],S[i])<0){
        tmp = S[i-1];
        S[i-1] = S[i];
        S[i] = tmp;
        bubblesort(S,n);
     }
     i++;
  }

  return;
}

/*****************************************************************************
@
@-----------------------------------------------------------------------------
@
@ char *symbol(bravais_TYP *G,matrix_TYP *F)
@
@ bravais_TYP *G: the group in question
@ matrix_TYP *F : a positive definite G-invariant form
@
@ Calculates the symbol of the FINITE INTEGRAL group in G.
@ F is assumed to be a positive definite G-invariant form.
@ The information needed from G are the records G->gen & G->gen_no,
@ and G->form & G->form_no. 
@
@ Sideefects: The matrices in G->gen and G->form might be checked via
@             Check_mat
@-----------------------------------------------------------------------------
@
*****************************************************************************/
char *symbol(bravais_TYP *G,matrix_TYP *F)
{
  int i,
      j,
      len=0,
      hom_no;      /* the number of homogenous constituents of this repr */

  constituent *hom;

  bravais_TYP *brav;  /* holds the bravais group of a homogenous
                         representation */

  char *result,
      **symb;

  /* firstly split the representation into homogenous parts */
  hom = homogenous(G,F,&hom_no);

  if (INFO_LEVEL & 4){
     for (i=0;i<hom_no;i++){
        put_bravais(hom[i].group,NULL,NULL);
     }
  }

  symb = (char **) calloc(hom_no , sizeof(char *));
  /* calculate the bravais group for each homogenous module,
     and search it in the list */
  for (i=0;i<hom_no;i++){
     brav = bravais_group(hom[i].group,FALSE);
     symb[i] = identify_hom(brav,i==(hom_no-1));
     free_bravais(brav);
     len += strlen(symb[i]);
  }

  /* swap the symbols in the right order */
  bubblesort(symb,hom_no);

  result = (char *) calloc((len + hom_no + 10),sizeof(char));
  sprintf(result,"%s",symb[0]);
  for (i=1;i<hom_no;i++){
     sprintf(result,"%s;%s",result,symb[i]);
  }

  for (i=0;i<hom_no;i++){
    free(symb[i]);
    free_bravais(hom[i].group);
    for (j=0;j<hom[i].dimc;j++)
       free_mat(hom[i].centralizer[j]);
    for (j=0;j<hom[i].dimcc;j++)
       free_mat(hom[i].ccentralizer[j]);

    if (hom[i].centralizer != NULL){
       free(hom[i].centralizer);
       hom[i].centralizer = NULL;
    }
    if (hom[i].ccentralizer != NULL){
       free(hom[i].ccentralizer);
       hom[i].ccentralizer = NULL;
    }

  }
  free(symb);
  free(hom);

  return result;

}

