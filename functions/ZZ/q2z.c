#include "typedef.h"
#include "voronoi.h"
#include "symm.h"
#include "autgrp.h"
#include "bravais.h"
#include "base.h"
#include "idem.h"
#include "longtools.h"
#include "reduction.h"
#include "ZZ_P.h"
#include "ZZ_zclass_P.h"
#include "ZZ_cen_fun_P.h"
#include "graph.h"
#include "datei.h"

int IDEM_NO;
extern int INFO_LEVEL;
ZZ_super_TYP **SUPER_info, *SUPER_INFO;


/* -------------------------------------------------------------------- */
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


   /* free konst
   if (flag == 1){
      for (i = 0; i < super->konst->k; i++){
         for (j = 0; j < super->konst->s[i]; j++){
            for (k = 0; k < super->konst->r; k++){
               free_mat(super->konst->Delta[i][j][k]);
            }
            free(super->konst->Delta[i][j]);
         }
         free(super->konst->Delta[i]);
      }
      free(super->konst->Delta);
      free(super->konst->s);
      free(super->konst);
   }
*/

/* ------------------------------------------------------------------ */
static void free_ZZ_super(ZZ_super_TYP *super)
{
   int i, j, k;
	
   ZZ_node_t *n, *t;


   n = super->tree->root;
   do {
      t = n->next;		
      ZZ_free_node (super->data, n);
   } while ((n = t) != NULL);
   free(super->tree);

   ZZ_free_data(super->data);
   free(super->data);

   free(super);
}




/* ------------------------------------------------------------------ */
void free_QtoZ(QtoZ_TYP *inz,
               int flag)
{
   int i, j, k;


   if (inz->zoogitter != NULL){
      for (i = 0; i < inz->anz; i++){
         for (j = 0; j < inz->entry[i][0].anz; j++){
            if (inz->zoogitter[i][j] != NULL)
               free_mat(inz->zoogitter[i][j]);
         }
         if (inz->zoogitter[i] != NULL)
            free(inz->zoogitter[i]);
      }
      free(inz->zoogitter);
   }
   for (i = 0; i < inz->anz; i++){
      if (flag == 1 && inz->gitter[i] != NULL){
         free_mat(inz->gitter[i]);
         free_mat(inz->tr_gitter[i]);
         free_mat(inz->inv_tr_gitter[i]);
      }
      for (j = 0; j < inz->anz + flag; j++){
         if (inz->entry[i][j].anz != 0){
            free(inz->entry[i][j].I);
            free(inz->entry[i][j].J);
            free(inz->entry[i][j].flag);
            if (flag == 0){
               for (k = 0; k < inz->entry[i][j].anz; k++){
                  free_mat(inz->entry[i][j].lattice[k]);
               }
            }
            free(inz->entry[i][j].lattice);
            if (flag == 0){
               for (k = 0; k < inz->entry[i][j].anz; k++){
                  free_mat(inz->entry[i][j].lsf[k]);
               }
               free(inz->entry[i][j].lsf);
            }
         }
      }
      free(inz->entry[i]);
   }
   free(inz->entry);

   if (flag == 1){
      free(inz->gitter);
      free(inz->inv_tr_gitter);
      free(inz->tr_gitter);
   }

   free(inz);
}



/*------------------------------------------------------------------------------- */
static int suche_mat(matrix_TYP *mat,
                     matrix_TYP **liste,
                     int anz)
{
   int i;

   for (i = 0; i < anz; i++){
      if (cmp_mat(mat, liste[i]) == 0)
         return(i);
   }
   return(-1);
}



/* ------------------------------------------------------------------------------- */
static matrix_TYP *hom_de_super(bravais_TYP *group,
                                matrix_TYP *id,
                                matrix_TYP **idem,
                                int idem_no)
{
   matrix_TYP *F,
              *tmp,
              *new_basis,
             **idem_spaces;

   int i, j, k,
       col, kgv,
       mult;



   F = rform(group->gen, group->gen_no, id, 101);
   idem_spaces = (matrix_TYP **)calloc(idem_no, sizeof(matrix_TYP *));

   for (i = 0; i < idem_no; i++){
      tmp = tr_pose(idem[i]);
      idem_spaces[i] = ggauss(tmp);
      idem_spaces[i] = good_initial_basis(idem_spaces[i], F);
      free_mat(tmp);
   }

   kgv = 1;
   for (i = 0; i < idem_no; i++){
      kgv = KGV(kgv, idem_spaces[i]->kgv);
   }

   for (i = 0; i < idem_no; i++){
      mult = kgv/idem_spaces[i]->kgv;
      if (mult != 1){
         for (j = 0; j < idem_spaces[i]->rows; j++){
            for (k = 0; k < idem_spaces[i]->cols; k++)
               idem_spaces[i]->array.SZ[j][k] *= mult;
         }
         idem_spaces[i]->kgv = kgv;
      }
   }

   new_basis = init_mat(group->dim, group->dim, "1");
   col = -1;
   for (i = 0; i < idem_no; i++){
      for (j = 0; j < idem_spaces[i]->rows; j++){
         col++;
         for (k = 0; k < group->dim; k++){
            new_basis->array.SZ[k][col] = idem_spaces[i]->array.SZ[j][k];
         }
      }
   }
   new_basis->kgv = kgv;
   Check_mat(new_basis);

   /* clean */
   for (i = 0; i < idem_no; i++){
      free_mat(idem_spaces[i]);
   }
   free(idem_spaces);
   free_mat(F);

   return(new_basis);
}



/* ------------------------------------------------------------------ */
static void get_normalizer(bravais_TYP *H)
{

    matrix_TYP *ID,
               *F,
               *PF,
              **N,
              **forms,
               *trbifo;

    bravais_TYP *Htr,
		  *HB;

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
	 HB = bravais_group(H,TRUE);
    V = normalizer(PF,HB,Htr,1949,&vno);
	 H->normal = HB->normal ; H->normal_no = HB->normal_no;
	 HB->normal = NULL; HB->normal_no = 0;
	 free_bravais(HB);

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



/* ------------------------------------------------------------------ */
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

/* changed by oliver:
   for (i=0;i<zentr_no;i++){
      free_mat(H->zentr[i]);
   }
   free(H->zentr);
   H->zentr = NULL;
*/
   H->zentr_no = zentr_no;

   return RES;
}



/* ------------------------------------------------------------------ */
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



/* ------------------------------------------------------------------ */
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
   mat_muleq(B, TR);

   return(TR);
}


static matrix_TYP **get_better_base(bravais_TYP *H,
                                    matrix_TYP *F,
                                    int number,
                                    matrix_TYP **SPACES)
{

   int i;

   matrix_TYP *id = NULL,
              **better,
               *tmp;


   if (F == NULL){
      id = init_mat(H->dim,H->dim,"i1");
      F = rform(H->gen,H->gen_no,id,101);
   }
   better = (matrix_TYP **)calloc(H->zentr_no, sizeof(matrix_TYP *));

   for (i = 1; i < H->zentr_no; i++){
      better[i] = better_base(H->zentr[i],F,number,SPACES);
   }

   if (id){
      free_mat(F);
      free_mat(id);
   }

   return(better);

}




static matrix_TYP **GL_n_Z(int dim,int *no)
{

   int i;

   matrix_TYP **RES;

   if (dim == 1){
      *no = 1;
      RES = (matrix_TYP **) malloc(1*sizeof(matrix_TYP *));
      RES[0] = init_mat(1,1,"i0");
      RES[0]->array.SZ[0][0] = -1;
   }
   else{
      *no = 4;
      RES = (matrix_TYP **) malloc(4*sizeof(matrix_TYP *));
      RES[0] = init_mat(dim,dim,"i1");
      RES[0]->array.SZ[0][0] = -1;
      RES[1] = init_mat(dim,dim,"i1");
      RES[1]->array.SZ[0][1] = 1;
      RES[2] = init_mat(dim,dim,"i1");
      RES[2]->array.SZ[0][1] = 1;
      RES[2]->array.SZ[1][0] = 1;
      RES[2]->array.SZ[0][0] = 0;
      RES[2]->array.SZ[1][1] = 0;
      RES[3] = init_mat(dim,dim,"i0");
      for (i=0;i<dim;i++)
         RES[3]->array.SZ[i][(i+1)%dim] = 1;
      for (i=0;i<4;i++)
         Check_mat(RES[i]);
   }

   return RES;

}

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
bravais_TYP **q2z(bravais_TYP *G,
                  int *number,
                  int ADFLAG,
                  QtoZ_TYP *INZ,
                  int quiet)
{
   bravais_TYP **GROUPS,
               **ADGROUPS,
                *H,
                *group,
                *ggg;

   matrix_TYP **IDEM,
              **IDEM_SPACES,
              **idem_spaces,
               *F,
               *tmp,
               *new_base,
               *id,
               *X,
/*               *test, */
              **conj_idem,
              **real_idem,
               *zoo_kon,
               *lattice,
               *zoolattice,
               *zoolattice_hnf,
               *zoolattice_inv,
               *zoo_inv,
               *elementar,
              **better,
              **trash,
               *zwischen,
               *ttt,
              **zentr_inv;

   char zzoptions[128],
        string[128];

   int i,
       j,
       k,
       l,
       m,
       f,
       kgv, mult,
       aa, bb,
       counter,
       ganzzahlig,
       help,
       idem_no,
       ad_no,
       col,
       dimc,
       dimcc,
       zen_nr, cen_nr, normal_nr,
       zoo, bahnnummer,
       samezoo, yeah,
      *liste;

   QtoZ_TYP **inzidenz;

   ZZ_node_t *node, *Knoten;

   ZZ_couple_t *laeufer;


   /* handle the case of the groups G=<I> and G=<-I> differently,
      because the first one couldn't be handled by this anyway,
      and the second one trivial is but computationaly hard in dim 5 & 6 */
   k = TRUE; l = TRUE;
   for (i=0;i<G->gen_no && k;i++){
      Check_mat(G->gen[i]);
      k = k && G->gen[i]->flags.Scalar;
      l = l && G->gen[i]->array.SZ[1][1] == 1;
   }
   if (l && k && G->order == 1)
      for (i = 0;i<100 ; l = l && G->divisors[i] == 0, i++);
   if ((k && G->dim > 4) ||
      (k && l )){
      GROUPS = (bravais_TYP **) malloc(2 * sizeof(bravais_TYP *));
      GROUPS[0] = copy_bravais(G);
      GROUPS[1] = (bravais_TYP * ) 1;
      if (GROUPS[0]->normal && GROUPS[0]->normal_no > 0){
         for (k=0;k<GROUPS[0]->normal_no;k++) free_mat(GROUPS[0]->normal[i]);
         free(GROUPS[0]->normal);
      }
      GROUPS[0]->normal = GL_n_Z(G->dim,&GROUPS[0]->normal_no);
      *number = 1;
      if (ADFLAG){
          GROUPS[0]->zentr = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP));
          GROUPS[0]->zentr[0] = init_mat(G->dim,G->dim,"1");
          GROUPS[0]->zentr_no = 1;
      }
      return GROUPS;
   }

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
   kgv = 1;
   for (i = 0; i < IDEM_NO; i++){
      kgv = KGV(kgv, IDEM_SPACES[i]->kgv);
   }

   for (i = 0; i < IDEM_NO; i++){
      mult = kgv/IDEM_SPACES[i]->kgv;
      if (mult != 1){
         for (j = 0; j < IDEM_SPACES[i]->rows; j++){
            for (k = 0; k < IDEM_SPACES[i]->cols; k++)
               IDEM_SPACES[i]->array.SZ[j][k] *= mult;
         }
         IDEM_SPACES[i]->kgv = kgv;
      }
   }

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

   new_base->kgv = kgv;
   Check_mat(new_base);
   tmp = mat_inv(new_base);
   i = G->zentr_no;
   G->zentr_no = 0;
   H = konj_bravais(G,tmp);
   for (m = 0; m < H->form_no; m++){
      H->form[m]->kgv = 1;
   }
   long_rein_formspace(H->form, H->form_no, 0);
   real_idem = &H->gen[H->gen_no - IDEM_NO];
   G->zentr_no = i;
   free_mat(tmp);

   /* make the ZZ options for the first call of ZZ */
   if (quiet)
      sprintf(zzoptions,"qtugZ");
   else
      sprintf(zzoptions,"tugZ");

   if (INFO_LEVEL & 4){
      put_mat(new_base,NULL,"new_base",2);
      put_bravais(H,NULL,"H");
      printf("zzoptions:  %s\n",zzoptions);
   }

   /* call the ZZ to get all almost decomposable groups in this Q-class */
   free_mat(F);
   F = rform(H->gen,H->gen_no- IDEM_NO,id,101);
   ZZ(H,F,NULL,NULL,zzoptions,NULL,0,-1);

   if (INFO_LEVEL & 4){
      put_bravais(H,NULL,"H");
   }

   /* apply a base reduction algorithm to the found new basises
      BASISES ARE CHANGED !!!!!!!!!!!!!!!!!!!!!!!!!!! */
   if (GRAPH){
      better = get_better_base(H,F,IDEM_NO,IDEM_SPACES);
      better[0] = init_mat(G->dim,G->dim,"1");
   }
   else{
      trash = get_better_base(H,F,IDEM_NO,IDEM_SPACES);
      for (j = 1; j < H->zentr_no; j++){
         free_mat(trash[j]);
      }
      free(trash);
   }

   /* get all almost decomposable groups */
   ad_no = H->zentr_no;
   ADGROUPS = almost(H);

   if (INFO_LEVEL & 4){
      fprintf(stderr,"=== Number of almost decomposable groups %d\n",ad_no);
   }

   /* make the ZZ options for the second call of ZZ */
   if (quiet)
      sprintf(zzoptions,"tuqzgp%d",IDEM_NO);
   else
      sprintf(zzoptions,"tuzgp%d",IDEM_NO);

   for (i = 0; i < IDEM_NO; i++){
      sprintf(string,"%s/%d",zzoptions,IDEM_SPACES[i]->rows);
      sprintf(zzoptions,"%s",string);
   }

   /* call the ZZ again to get all sublattices which fullfill the projection
      property (for each group in  ADGROUPS) */
   idem_no = IDEM_NO;
   IDEM_NO = 0;

   if (GRAPH){
      inzidenz = (QtoZ_TYP **)calloc(ad_no, sizeof(QtoZ_TYP *));
      SUPER_info = (ZZ_super_TYP **)calloc(ad_no, sizeof(ZZ_super_TYP *));
      for (i = 0; i < ad_no; i++){
         inzidenz[i] = (QtoZ_TYP *)calloc(1, sizeof(QtoZ_TYP));
         inzidenz[i]->gitter = (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
         inzidenz[i]->inv_tr_gitter = (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
         inzidenz[i]->tr_gitter = (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
         inzidenz[i]->zoogitter = (matrix_TYP ***)calloc(1024, sizeof(matrix_TYP **));
         inzidenz[i]->entry = (QtoZ_entry_TYP **)calloc(1024, sizeof(QtoZ_entry_TYP *));
         ZZ(ADGROUPS[i],F,NULL,inzidenz[i],zzoptions,NULL,i,(i == 0));
      }
   }
   else{
      for (i = 0; i < ad_no; i++){
         ZZ(ADGROUPS[i],F,NULL,NULL,zzoptions,NULL,i,0);

	 /* apply a base reduction algorithm to the found new basises */
         trash = get_better_base(ADGROUPS[i],NULL,1,ADGROUPS[i]->gen);
         for (j = 1; j < ADGROUPS[i]->zentr_no; j++){
            free_mat(trash[j]);
         }
         free(trash);
      }
   }

   /* apply a base reduction algorithm to the found new basises */
/*
   for (i = 0; i < ad_no; i++){
      trash = get_better_base(ADGROUPS[i],NULL,1,ADGROUPS[i]->gen);
      for (j = 1; j < inzidenz[i]->anz; j++){
         free_mat(trash[j]);
      }
      free(trash);
   }
*/

   if (GRAPH){
      zentr_inv = (matrix_TYP **)calloc(ad_no, sizeof(matrix_TYP *));
      for (i = 0; i < ad_no; i++){
         zentr_inv[i] = mat_inv(H->zentr[i]);
         INZ->anz += inzidenz[i]->anz;
      }
      INZ->entry = (QtoZ_entry_TYP **)calloc(INZ->anz, sizeof(QtoZ_entry_TYP *));
      for (i = 0; i < INZ->anz; i++){
         INZ->entry[i] = (QtoZ_entry_TYP *)calloc(INZ->anz, sizeof(QtoZ_entry_TYP));
      }
      liste = (int *)calloc(ad_no + 1, sizeof(int));
      liste[0] = counter = 0;
      for (i = 0; i < ad_no; i++){
         for (j = 0; j < inzidenz[i]->anz; j++){
            for (l = 0 ; l < inzidenz[i]->anz; l++){
               help = INZ->entry[j + counter][l + counter].anz = inzidenz[i]->entry[j][l + 1].anz;
               if (help != 0){
                  INZ->entry[j + counter][l + counter].I = (int *)calloc(1024, sizeof(int));
                  INZ->entry[j + counter][l + counter].J = (int *)calloc(1024, sizeof(int));
                  INZ->entry[j + counter][l + counter].flag = (int *)calloc(1024, sizeof(int));
                  INZ->entry[j + counter][l + counter].lattice =
                                         (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
                  for (m = 0; m < help; m++){
                     INZ->entry[j + counter][l + counter].I[m] = inzidenz[i]->entry[j][l + 1].I[m];
                     INZ->entry[j + counter][l + counter].J[m] = inzidenz[i]->entry[j][l + 1].J[m];
                     INZ->entry[j + counter][l + counter].flag[m] =
                                              inzidenz[i]->entry[j][l + 1].flag[m];
                     INZ->entry[j + counter][l + counter].lattice[m] =
                                                           inzidenz[i]->entry[j][l + 1].lattice[m];
                     Check_mat(INZ->entry[j + counter][l + counter].lattice[m]);
                  }
               }
            }
         }
         counter += inzidenz[i]->anz;
         liste[i+1] = counter;
      }

      /* handle the lattices, which might be in another zoo */
      for (i = 0; i < ad_no; i++){
         zen_nr = ADGROUPS[i]->zentr_no;
         normal_nr = ADGROUPS[i]->normal_no;
         cen_nr = ADGROUPS[i]->cen_no;
         ADGROUPS[i]->zentr_no = 0;
         ADGROUPS[i]->normal_no = 0;
         ADGROUPS[i]->cen_no = 0;
         for (j = 0; j < inzidenz[i]->anz; j++){
            for (l = 0; l < inzidenz[i]->entry[j][0].anz; l++){
               samezoo = 0;

               /* calculate group */
               zoo_inv = mat_inv(inzidenz[i]->zoogitter[j][l]);
               ggg = konj_bravais(ADGROUPS[i], zoo_inv);
               for (m = 0; m < ggg->form_no; m++){
                  ggg->form[m]->kgv = 1;
               }
               long_rein_formspace(ggg->form, ggg->form_no, 0);

               /* conjugate idempotents (the idempotents are equal for all
                  homogeneously decomposable lattices) */
               conj_idem = (matrix_TYP **)calloc(idem_no, sizeof(matrix_TYP *));
               ganzzahlig = 1;
               for (m = 0; m < idem_no; m++){
                  conj_idem[m] = mat_kon(zoo_inv ,real_idem[m], inzidenz[i]->zoogitter[j][l]);
                  Check_mat(conj_idem[m]);
                  if (conj_idem[m]->kgv != 1){
                     ganzzahlig = 0;
                  }
               }

               if (ganzzahlig == 1){
                  /* it's an homogeneously decomposable lattice */
                  group = ggg;
               }
               else{
                  /* calculate the minimal homogeneously decomposable superlattice */
                  lattice = hom_de_super(ggg, id, conj_idem, idem_no);
                  tmp = mat_inv(lattice);
                  group = konj_bravais(ggg, tmp);
	          for (f = 0; f < group->form_no; f++)
		     group->form[f]->kgv = 1;
	          long_rein_formspace(group->form, group->form_no, 1);
               }

               /* which one of the old homogeneously decomposable lattices is it */
               zoo_kon = special_deal_with_zclass(SUPER_INFO->tree, group, &zoo);

               if (ganzzahlig != 1){
                  if (zoo != i){
                     /* another zoo */
                     samezoo = 0;
                  }
                  else{
                     /* same zoo */
                     samezoo = 1;
                  }
                  Knoten = SUPER_INFO->tree->root;
                  for (f = 0; f < zoo; f++){
                     Knoten = Knoten->next;
                  }
                  zwischen = mat_kon(zoo_kon, Knoten->Q, better[zoo]);
                  ttt = mat_inv(zwischen);
                  zoolattice = mat_mul(ttt, tmp);
                  free_mat(ttt);
                  Knoten = NULL;
               }

               /* free */
               if (ganzzahlig == 1){
                  group = NULL;
               }
               else{
                  free_bravais(group);
               }
               for (m = 0; m < idem_no; m++){
                  free_mat(conj_idem[m]);
               }
               free(conj_idem);
               free_bravais(ggg);

               /* which lattice in the zoo? */
               aa = liste[i] + j;

               if (ganzzahlig == 1){
                  /* it is the homogeneously decomposable lattice */
                  bb = liste[zoo];
                  bahnnummer = 0;
               }
               else{
                  /* search in the tree */
                  zoolattice_hnf = copy_mat(zoolattice);
                  long_col_hnf(zoolattice_hnf);
                  node = SUPER_info[zoo]->tree->root;
                  yeah = 0;
                  while (node != NULL){
                     /* could be made better: compare el_div first */
                     if (in_bahn(zoolattice_hnf, node, &bahnnummer) == 1){
                        yeah = 1;
                        break;
                     }
                     node = node->next;
                  }
                  if (yeah == 0){
                     fprintf(stderr,"ERROR 1 in q2z!\n");
                     exit(3);
                  }
  	          laeufer = node->child;
      	          for (m = 0; m < node->N_no_orbits - bahnnummer; m++){
	             laeufer = laeufer->elder;
  	             if (laeufer == NULL){
   	                fprintf(stderr,"ERROR 2 in q2z!\n");
 	                exit(2);
 	             }
	          }
  	          bahnnummer = suche_mat(laeufer->he->U, inzidenz[zoo]->gitter, inzidenz[zoo]->anz);
                  if (bahnnummer == -1){
	             fprintf(stderr,"ERROR 3 in q2z!\n");
 	             exit(3);
	          }
	          bb = liste[zoo] + bahnnummer;
	          node = NULL;

                  free_mat(zoolattice_hnf);
               }

               /* write information in INZ */
               help = INZ->entry[aa][bb].anz;
               if (help == 0){
                  INZ->entry[aa][bb].I = (int *)calloc(1024, sizeof(int));
                  INZ->entry[aa][bb].J = (int *)calloc(1024, sizeof(int));
                  INZ->entry[aa][bb].flag = (int *)calloc(1024, sizeof(int));
                  INZ->entry[aa][bb].lattice = (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
               }
               INZ->entry[aa][bb].I[help] = inzidenz[i]->entry[j][0].I[l];
               INZ->entry[aa][bb].J[help] = inzidenz[i]->entry[j][0].J[l];
               INZ->entry[aa][bb].flag[help] = inzidenz[i]->entry[j][0].flag[l];

               if (ganzzahlig == 1){
                  Knoten = SUPER_INFO->tree->root;
                  for (f = 0; f < zoo; f++){
                     Knoten = Knoten->next;
                  }
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], inzidenz[i]->zoogitter[j][l]);
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], zoo_kon);
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], Knoten->Q);
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], better[zoo]);
               }
               else{
                  zoolattice_inv = mat_inv(zoolattice);
                  X = konjugierende(zoolattice_inv, SUPER_info[zoo]->tree->root->col_group,
                                    laeufer->he);
                  free_mat(zoolattice_inv);
                  if (X == NULL){
                     fprintf(stderr, "ERROR 4 in q2z!\n");
                     exit(9);
                  }

                  /*
                  if (samezoo != 1){
                     mat_muleq(inzidenz[i]->entry[j][0].lattice[l], zentr_inv[j]);
                     mat_muleq(inzidenz[i]->entry[j][0].lattice[l], H->zentr[zoo]);
                  }
                  */
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], inzidenz[i]->zoogitter[j][l]);
                  /*
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], lattice);
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], zwischen);
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], zoolattice);
                  */
                  mat_muleq(inzidenz[i]->entry[j][0].lattice[l], X);
                  free_mat(X);
                  free_mat(zoolattice);
                  free_mat(lattice);
                  free_mat(tmp);
                  free_mat(zwischen);
               }
               INZ->entry[aa][bb].lattice[help] = inzidenz[i]->entry[j][0].lattice[l];
               Check_mat(INZ->entry[aa][bb].lattice[help]);
               INZ->entry[aa][bb].anz++;
               if (zoo_inv != NULL)
                  free_mat(zoo_inv);
               free_mat(zoo_kon);
            }
         }
         ADGROUPS[i]->zentr_no = zen_nr;
         ADGROUPS[i]->normal_no = normal_nr;
         ADGROUPS[i]->cen_no = cen_nr;
      }

      /* print the matrix of incidences for the Z-classes */
      /*
      test = init_mat(INZ->anz, INZ->anz, "");
      for (i = 0; i < INZ->anz; i++){
         for (j = 0; j < INZ->anz; j++){
            test->array.SZ[i][j] = INZ->entry[i][j].anz;
         }
      }
      put_mat(test,0,"graph for the arithmetic classes",0);
      free_mat(test);
      */

      /* calculate standard form for the lattices */
      for (i = 0; i < INZ->anz; i++){
         for (j = 0; j < INZ->anz; j++){
            if (INZ->entry[i][j].anz != 0){
               INZ->entry[i][j].lsf = (matrix_TYP **)calloc(INZ->entry[i][j].anz,
                                                            sizeof(matrix_TYP *));
               for (k = 0; k < INZ->entry[i][j].anz; k++){
                  INZ->entry[i][j].lsf[k] = copy_mat(INZ->entry[i][j].lattice[k]);
                  long_col_hnf(INZ->entry[i][j].lsf[k]);
               }
            }
         }
      }
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

   if (GRAPH){
      /* free better, inzidenz, SUPER_info and SUPER_INFO (oliver: 9.8.00) */
      free_ZZ_super(SUPER_INFO);
      for (i = 0; i < ad_no; i++){
         free_mat(zentr_inv[i]);
         free_ZZ_super(SUPER_info[i]);
         free_QtoZ(inzidenz[i], 1);
         free_mat(better[i]);
      }
      free(zentr_inv);
      free(liste);
      free(SUPER_info);
      free(inzidenz);
      free(better);
      ZZ_get_data (NULL, NULL, NULL, NULL, NULL, NULL, -3);
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
