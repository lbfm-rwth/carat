#include "ZZ.h"
#include "typedef.h"
#include "voronoi.h"
#include "symm.h"
#include "autgrp.h"
#include "bravais.h"
#include "base.h"
#include "idem.h"
#include "longtools.h"
#include "reduction.h"

int IDEM_NO;
extern int INFO_LEVEL;

static void get_normalizer(bravais_TYP *H)
{

    matrix_TYP *ID,
               *F,
               *PF,
              **N,
              **forms,
               *trbifo;

    bravais_TYP *Htr;

    int i,
        j,
        vno;

    voronoi_TYP **V;

    Htr = tr_bravais(H,1,FALSE);
    ID = init_mat(H->dim,H->dim,"1");
    F = rform(H->gen,H->gen_no,ID,101);
    trbifo = trace_bifo(H->form,Htr->form,H->form_no);
    PF = first_perfect(F,H,Htr->form,trbifo,&i);

    /* calculate the normalizer of the bravais group */
    V = normalizer(PF,H,Htr,1949,&vno);

    /* calculate the bravais group itself */
    free_bravais(Htr);
    forms = (matrix_TYP **) malloc( (1+H->form_no) * sizeof(matrix_TYP *));
    for (i=0;i<H->form_no;i++)
       forms[i+1] = H->form[i];
    forms[0] = PF;     /* a H-perfect form is especially good */


    /* replaced the following two statements by a better one
    tilman 06.05.98:
    SV = short_vectors(PF,max_diagonal_entry(PF),0,0,0,&i);
    Htr = autgrp(forms,H->form_no+1,SV,NULL,0,NULL); */
    Htr = pr_aut(forms,H->form_no+1,NULL,0,NULL);

    Htr->normal = H->normal;
    Htr->normal_no = H->normal_no;
    H->normal = NULL;
    H->normal_no = 0;

    H->normal = normalizer_in_N(H,Htr,&H->normal_no,FALSE);

    /* clean up */
    free_mat(ID);
    free_mat(F);
    free_mat(PF);
    free_mat(trbifo);
    free_bravais(Htr);
    for (i=0;i<vno;i++){
       clear_voronoi(V[i]);
       free(V[i]);
    }
    free(V);
    /* free_mat(SV); */
    free(forms);

    return;
}

static bravais_TYP **almost(bravais_TYP *H)
{
   int i,
       j,
       zentr_no,
       normal_no,
       cen_no;

   matrix_TYP *tmp;

   bravais_TYP **RES;

   RES = (bravais_TYP **) malloc((1+H->zentr_no) * sizeof(bravais_TYP *));

   H->gen_no -= IDEM_NO;
   normal_no = H->normal_no;
   H->normal_no = 0;
   cen_no = H->cen_no;
   H->cen_no = 0;
   zentr_no = H->zentr_no;
   H->zentr_no = 0;
   for (i=0;i<zentr_no;i++){
      if (i==0){
         RES[i] = copy_bravais(H);
      }
      else{
         tmp = mat_inv(H->zentr[i]);
         RES[i] = konj_bravais(H,tmp);
         free_mat(tmp);
      }
      for (j=0;j<H->form_no;j++){
         RES[i]->form[j]->kgv = 1;
      }
      long_rein_formspace(RES[i]->form,RES[i]->form_no,1);
      get_normalizer(RES[i]);
   }

   H->gen_no += IDEM_NO;
   H->cen_no = cen_no;
   H->normal_no = normal_no;

   for (i=0;i<zentr_no;i++){
      free_mat(H->zentr[i]);
   }
   free(H->zentr);
   H->zentr = NULL;

   return RES;
}

static bravais_TYP **get_groups(bravais_TYP **ADGROUPS,int ad_no,int *number)
{
   int i,
       j,
       k,
       normal_no,
       centr_no;

   bravais_TYP **GROUPS,
                *G1,
                *G2;

   matrix_TYP  *tmp,
              **centerings,
              **normal;

   G1 = init_bravais(ADGROUPS[0]->dim);

   /* look for all the centerings and so on */
   number[0] = 0;
   for (i=0;i<ad_no;i++){
     number[0] += ADGROUPS[i]->zentr_no;
   }
   GROUPS = (bravais_TYP **) malloc(number[0] * sizeof(bravais_TYP *));

   k=0;
   for (i=0;i<ad_no;i++){
      centerings = ADGROUPS[i]->zentr;
      centr_no = ADGROUPS[i]->zentr_no;
      ADGROUPS[i]->zentr = NULL;
      ADGROUPS[i]->zentr_no = 0;
      normal = ADGROUPS[i]->normal;
      ADGROUPS[i]->normal = NULL;
      normal_no = ADGROUPS[i]->normal_no;
      ADGROUPS[i]->normal_no = 0;
      G1->gen = normal;
      G1->gen_no = normal_no;
      for (j=0;j<centr_no;j++){
         tmp = mat_inv(centerings[j]);
         if (j==0){
            GROUPS[k] = konj_bravais(ADGROUPS[i],tmp);
            GROUPS[k]->normal = normal;
            GROUPS[k]->normal_no = normal_no;
         }
         else{
            G2 = gittstab(G1,centerings[j]);
            ADGROUPS[i]->normal = G2->gen;
            ADGROUPS[i]->normal_no = G2->gen_no;
            GROUPS[k] = konj_bravais(ADGROUPS[i],tmp);
            long_rein_formspace(GROUPS[k]->form,GROUPS[k]->form_no,1);
            ADGROUPS[i]->normal = NULL;
            ADGROUPS[i]->normal_no = 0;
            free_bravais(G2);
         }
         free_mat(tmp);
         free_mat(centerings[j]);
         k++;
      }
      free(centerings);
   }

   free(G1);

   return GROUPS;
}

/* this function was for test purposes only
static matrix_TYP *better_base2(matrix_TYP *B,
                                int number,
                                matrix_TYP **SP)
{

   int i,
       j,
       k,
       offset = 0;

   matrix_TYP *NEW,
              *C;

   put_mat(B,NULL,"B",0);
   NEW = init_mat(B->rows,B->cols,"i");

   for (i=0;i<number;i++){
      C = init_mat(SP[i]->rows,B->rows,"i");
      for (j=0;j<SP[i]->rows;j++){
         for (k=0;k<B->rows;k++){
            C->array.SZ[j][k] = B->array.SZ[k][j+offset];
         }
      }

      long_row_basis(C,TRUE);

      for (j=0;j<SP[i]->rows;j++){
         for (k=0;k<B->rows;k++){
            NEW->array.SZ[k][j+offset] = C->array.SZ[j][k];
         }
      }
      free_mat(C);

      offset += SP[i]->rows;
   }

   free_mat(B);

   return NEW;
}  */

static matrix_TYP *better_base(matrix_TYP *B,
                               matrix_TYP *F,
                               int number,
                               matrix_TYP **SP)
{

   int i,
       j,
       k,
       offset = 0;


   matrix_TYP *NEW,
              *TR,
              *C,
              *D;

   /* put_mat(B,NULL,"B",0); */

   TR = init_mat(B->cols,B->cols,"i");
   for (i=0;i<number;i++){
      C = init_mat(SP[i]->rows,B->rows,"i");
      for (j=0;j<B->rows;j++){
         for (k=0;k<SP[i]->rows;k++){
            C->array.SZ[k][j] = B->array.SZ[j][k+offset];
         }
      }

      /* calculate the form on this subspace */
      NEW = scal_pr(C,F,TRUE);
      free_mat(C);

      /* reduce the form */
      C = init_mat(NEW->rows,NEW->rows,"i1");
      D = pair_red(NEW,C);
      free_mat(D);

      /* store the transformation */
      for (j=0;j<C->rows;j++){
         for (k=0;k<C->cols;k++){
            TR->array.SZ[offset+k][offset+j] = C->array.SZ[j][k];
         }
      }

      free_mat(C);
      free_mat(NEW);
      offset += SP[i]->rows;
   }

   /* put_mat(TR,0,"TR",0); */
   NEW = mat_mul(B,TR);
   free_mat(B);
   free_mat(TR);

   return NEW;
}


static void get_better_base(bravais_TYP *H,
                            matrix_TYP *F,
                            int number,
                            matrix_TYP **SPACES)
{

   int i;

   matrix_TYP *id = NULL;

   if (F == NULL){
      id = init_mat(H->dim,H->dim,"i1");
      F = rform(H->gen,H->gen_no,id,101);
   }

   for (i=1;i<H->zentr_no;i++){
      H->zentr[i] = better_base(H->zentr[i],F,number,SPACES);
   }

   if (id){
      free_mat(F);
      free_mat(id);
   }

   return;

}

static matrix_TYP *good_initial_basis(matrix_TYP *B,
                                      matrix_TYP *F)
{

   int i;

   matrix_TYP *NEW,
              *C,
              *D;

   NEW = scal_pr(B,F,TRUE);

   /* reduce */
   C = init_mat(NEW->rows,NEW->rows,"i1");
   D = pair_red(NEW,C);

   free_mat(D);
   free_mat(NEW);

   NEW = mat_mul(C,B);

   free_mat(C);
   free_mat(B);

   return NEW;
}

bravais_TYP **q2z(bravais_TYP *G,int *number,int ADFLAG)
{
   bravais_TYP **GROUPS,
               **ADGROUPS,
                *H;

   matrix_TYP **IDEM,
              **IDEM_SPACES,
               *F,
               *tmp,
               *new_base,
               *id;

   char zzoptions[128],
        string[128];

   int i,
       j,
       k,
       idem_no,
       ad_no,
       col,
       dimc,
       dimcc;

   /* avoid silly mistakes */
   if (G->form == NULL ||
      G->form_no == 0 ||
      G->order == 0 ){
      fprintf(stderr,"You didn't specify G in an accurate way\n");
      exit(3);
   }

   id = init_mat(G->dim,G->dim,"1");

   /* we need a G-invariant, positive definite form for various reasons */
   F = rform(G->gen,G->gen_no,id,101);

   /* get the idempotents of the group */
   IDEM = idempotente(G->gen,G->gen_no,F,&IDEM_NO,&dimc,&dimcc,NULL);
   G->gen = (matrix_TYP **) realloc(G->gen,(G->gen_no+IDEM_NO)
                                    * sizeof(matrix_TYP*));
   IDEM_SPACES = (matrix_TYP **) malloc(IDEM_NO * sizeof(matrix_TYP *));
   for (i=0;i<IDEM_NO;i++){
      G->gen[i+G->gen_no] = IDEM[i];
      tmp = tr_pose(IDEM[i]);
      IDEM_SPACES[i] = long_rein_mat(tmp);
      IDEM_SPACES[i] = good_initial_basis(IDEM_SPACES[i],F);
      free_mat(tmp);
   }
   G->gen_no += IDEM_NO;

   /* transform G to H which is an almost decomposable group */
   new_base = init_mat(G->dim,G->dim,"1");
   col = -1;
   for (i=0;i<IDEM_NO;i++){
      for (j=0;j<IDEM_SPACES[i]->rows;j++){
         col++;
         for (k=0;k<G->dim;k++){
            new_base->array.SZ[k][col] = IDEM_SPACES[i]->array.SZ[j][k];
         }
      }
   }

   Check_mat(new_base);
   tmp = mat_inv(new_base);
   i = G->zentr_no;
   G->zentr_no = 0;
   H = konj_bravais(G,tmp);
   G->zentr_no = i;
   free_mat(tmp);

   /* make the ZZ options for the first call of ZZ */
   sprintf(zzoptions,"tugZ");

   if (INFO_LEVEL & 4){
      put_mat(new_base,NULL,"new_base",2);
      put_bravais(H,NULL,"H");
      printf("zzoptions:  %s\n",zzoptions);
   }

   /* call the ZZ to get all almost decomposable groups in this Q-class */
   free_mat(F);
   F = rform(H->gen,H->gen_no- IDEM_NO,id,101);
   ZZ(H,F,NULL,zzoptions,NULL);

   if (INFO_LEVEL & 4){
      put_bravais(H,NULL,"H");
   }

   /* apply a base reduction algorithm to the found new basises */
   get_better_base(H,F,IDEM_NO,IDEM_SPACES);

   /* get all almost decomposable groups */
   ad_no = H->zentr_no;
   ADGROUPS = almost(H);

   if (INFO_LEVEL & 4){
   }
      fprintf(stderr,"=== Number of almost decomposable groups %d\n",ad_no);

   /* make the ZZ options for the second call of ZZ */
   sprintf(zzoptions,"tuzgp%d",IDEM_NO);
   for (i=0;i<IDEM_NO;i++){
      sprintf(string,"%s/%d",zzoptions,IDEM_SPACES[i]->rows);
      sprintf(zzoptions,"%s",string);
   }

   /* call the ZZ again to get all sublattices which fullfill the projection
      property (for each group in  ADGROUPS) */
   idem_no = IDEM_NO;
   IDEM_NO = 0;
   for (i=0;i<ad_no;i++){
      ZZ(ADGROUPS[i],F,NULL,zzoptions,NULL);
      /* apply a base reduction algorithm to the found new basises */
      get_better_base(ADGROUPS[i],NULL,1,ADGROUPS[i]->gen);

   }

   if (!ADFLAG){
      GROUPS = get_groups(ADGROUPS,ad_no,number);
   }

   /* clean up */
   if (!ADFLAG){
      for (i=0;i<ad_no;i++){
         free_bravais(ADGROUPS[i]);
      }
      free(ADGROUPS);
   }
   free_bravais(H);
   free_mat(id);
   free_mat(new_base);
   free_mat(F);
   for (i=0;i<idem_no+dimc+dimcc;i++){
      free_mat(IDEM[i]);
   }
   free(IDEM);
   for (i=0;i<idem_no;i++){
      free_mat(IDEM_SPACES[i]);
   }
   free(IDEM_SPACES);
   G->gen_no -= idem_no;

   if (ADFLAG){
      number[0] = ad_no;
      return ADGROUPS;
   }
   else{
      return GROUPS;
   }
}
