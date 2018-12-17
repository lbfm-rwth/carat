#include"typedef.h"
#include"matrix.h"
#include"idem.h"
#include"gmp.h"
#include"symm.h"
#include"bravais.h"
#include"longtools.h"
#include "getput.h"
#include "sort.h"

extern int INFO_LEVEL;

/**************************************************************************
@
@--------------------------------------------------------------------------
@--------------------------------------------------------------------------
@ FILE: centr.c
@--------------------------------------------------------------------------
@--------------------------------------------------------------------------
@
***************************************************************************/

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ matrix_TYP **solve_endo(matrix_TYP **A,matrix_TYP **B,int anz,int *dim)
@
@ Let A[0],..,A[anz-1] be mxm matrices, and B[0],..,B[anz-1] nxn matrices.
@ The function returns a Q-basis of the space of mxn matrices with
@ A[i] * X = X * B[i] for all i. The dimension is returned via dim[0].
@
@ SIDEEFFECTS: The matrices in A,B are checked with Check_mat and
@              converted with rat2kgv.
@--------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP **solve_endo(matrix_TYP **A,matrix_TYP **B,int anz,int *dim)
{
  int i,
      j,
      k,
      l,
      m,
      act_row,
      Arows=A[0]->rows,
      Brows=B[0]->rows;

  matrix_TYP *GLS,      /* the linar system of equation to be solved */
            **sols,     /* the solutions */
            **rows;     /* the solution as rows */

  if (INFO_LEVEL & 4){
     fprintf(stderr,"entered solve_endo\n");
  }

  /* check some trivialities to avoid trouble */
  for (i=0;i<anz;i++){
     if ((A[i]->cols != Arows) || (A[i]->cols != Arows) ||
         (B[i]->cols != Brows) || (B[i]->rows != Brows)){
         fprintf(stderr,"error in solve_endo\n");
         exit(3);
     }
     rat2kgv(A[i]);
     Check_mat(A[i]);
     rat2kgv(B[i]);
     Check_mat(B[i]);

     /* it shouldn't matter whether or not the matrices are rational
     if (A[i]->kgv != 1){
        fprintf(stderr,"error in solve_endo\n");
        exit(3);
     }
     if (B[i]->kgv != 1){
        fprintf(stderr,"error in solve_endo\n");
        exit(3);
     } */
  }

  /* set the linear system of equations, no big entries will occur */
  dim[0] = 0;
  GLS = init_mat(anz*Arows*Brows,Arows*Brows,"");

  for (i=0;i<anz;i++){
     for (j=0;j<Arows;j++){
        for (k=0;k<Brows;k++){

           /* set the (i*Arows*Brows + j*Brows + k)-th row of the
              system of equations */
           act_row = i*Arows*Brows + j*Brows + k;
           for (l=0;l<Arows;l++){
              for (m=0;m<Brows;m++){
                 if (m==k && A[i]->array.SZ[j][l] !=0)
                    GLS->array.SZ[act_row][l*Brows+m] = A[i]->array.SZ[j][l];
                 if (l==j && B[i]->array.SZ[m][k] !=0)
                    GLS->array.SZ[act_row][l*Brows+m] -= B[i]->array.SZ[m][k];
              }
           }
        }
     }
  }

  Check_mat(GLS);

  if (INFO_LEVEL & 4){
     put_mat(GLS,NULL,"GLS",2);
  }

  /* we got the relevant system of equations, solve it, and use
     multiple precision integer to avoid overflow */
  rows = long_solve_mat(GLS,NULL);

  /* we won't need the system of equations anymore */
  free_mat(GLS);

  if (rows[0] != NULL){
     fprintf(stderr,"error in solve_endo, rows[0] should be NULL\n");
     free_mat(rows[0]);
  }

  if (rows[1] != NULL){ 
     /* there is a solution */
     dim[0] = rows[1]->cols;
     sols = (matrix_TYP **) malloc(dim[0] * sizeof(matrix_TYP *));

     for (i=0;i<rows[1]->cols;i++){
        sols[i] = init_mat(Arows,Brows,"");
        for (l=0;l<Arows;l++)
           for (m=0;m<Brows;m++)
              sols[i]->array.SZ[l][m] = rows[1]->array.SZ[l*Brows+m][i];
     }
     free_mat(rows[1]);
  }
  else{
     sols = NULL;
  }

  free(rows);

  return sols;
}

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ matrix_TYP **calc_ccentralizer(matrix_TYP **centr,int dim_centr,
@                                matrix_TYP **gen,int gen_no,int *dim_cc)
@
@
@--------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP **calc_ccentralizer(matrix_TYP **centr,int dim_centr,
                               matrix_TYP **gen,int gen_no,int *dim_cc)
{
   int anz = dim_centr + gen_no;

   matrix_TYP **tmp,
              **erg;

   tmp = (matrix_TYP **) malloc(anz * sizeof(matrix_TYP *));

   /* the first matrices of tmp will be those of centr */
   memcpy(tmp,centr,dim_centr * sizeof(matrix_TYP*));

   /* the next will be those of gen */
   memcpy(tmp+dim_centr,gen,gen_no * sizeof(matrix_TYP*));

   erg = solve_endo(tmp,tmp,anz,dim_cc);

   free(tmp);

   return erg;
   
}


int is_zero(matrix_TYP *pol,int x)
{
   MP_INT sum,
          pot,
          tmp;

   int i;

   mpz_init_set_ui(&sum,0);
   mpz_init_set_ui(&pot,1);
   mpz_init(&tmp);

   for (i=0;i<pol->cols;i++){
      mpz_set_si(&tmp,pol->array.SZ[0][i]);
      mpz_mul(&tmp,&pot,&tmp);
      mpz_add(&sum,&sum,&tmp);
      mpz_set_si(&tmp,x);
      mpz_mul(&pot,&pot,&tmp);
   }

   if (mpz_cmp_si(&sum,0) == 0)
      i = TRUE;
   else
      i = FALSE;

   mpz_clear(&sum);
   mpz_clear(&pot);
   mpz_clear(&tmp);

   return i;

}

void div_by_root(matrix_TYP *pol,int x)
{
   int i,
      *l;

   l = (int *) malloc((pol->cols-1)*sizeof(int));

   l[pol->cols-2] = pol->array.SZ[0][pol->cols-1];

   for (i=pol->cols-3;i>=0 ;i--){
      l[i] =  pol->array.SZ[0][i+1] + l[i+1] * x;
   }

   pol->cols--;
   free(pol->array.SZ[0]);
   pol->array.SZ[0] = l;

}

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ matrix_TYP *zeros(matrix_TYP *minpol)
@
@---------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP *zeros(matrix_TYP *minpol)
{
  int i,
      j,
      found = 0;

  matrix_TYP *erg,
             *tmp;

  erg = init_mat(1,minpol->cols,"");
  tmp = copy_mat(minpol);

  /* split of all root which are zero */
  while(tmp->array.SZ[0][0] == 0){
     div_by_root(tmp,0);
     erg->array.SZ[0][found] = 0;
     found++;
  }

  /* and now all others */
  for (i=1;i<=abs(tmp->array.SZ[0][0]);i++){
     if (tmp->array.SZ[0][0] % i == 0){
        for (j= -1;j<=1;j+=2){      /* describes the sign */
           if (is_zero(tmp,j*i)){
              div_by_root(tmp,i*j);
              erg->array.SZ[0][found] = i * j;
              found++;
              i=0;
              j=3;
           }
        }
     }
  }

  free_mat(tmp);
  real_mat(erg,1,found);

  return erg;

} /* zeros(...) */

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ int pos(matrix_TYP **list,int no,matrix_TYP *x)
@
@--------------------------------------------------------------------------
@
***************************************************************************/
int pos(matrix_TYP **list,int no,matrix_TYP *x)
{
  int i;

  for (i=0;i<no;i++)
     if (mat_comp(list[i],x) == 0) return i;

  return -1;
} /* pos (....) */

static matrix_TYP *matrix_in_pol(matrix_TYP *pol,matrix_TYP *x)
{
   int i;

   matrix_TYP *erg,
              *pot;

   if (x->cols != x->rows && x->flags.Integral == FALSE){
      fprintf(stderr,"exit in matrix_in_pol\n");
      exit(3);
   }

   pot = init_mat(x->cols,x->cols,"1");
   erg = init_mat(x->cols,x->cols,"");

   for (i=0;i<pol->cols;i++){
      imat_addeq(erg,pot,1,pol->array.SZ[0][i]);
      mat_muleq(erg,x);
   }

   free_mat(pot);

   return erg;

} /* matrix_in_pol(...) */

int min_trace(matrix_TYP **A,int no)
{

  int i,
      min,
      p,
      trac;

  min = trace(A[0]);
  p = 0;
  if (min % A[0]->kgv != 0){
     fprintf(stderr,"error in min_trace\n");
     exit(3);
  }
  else{
     min = abs(min/A[0]->kgv);
  }

  for (i=1;i<no;i++){
     trac = trace(A[i]);
     if (trac % A[i]->kgv != 0){
        fprintf(stderr,"error in min_trace\n");
        exit(3);
     }
     else{
        trac = abs(trac/A[i]->kgv);
     }
     if (trac < min){
        min = trac;
        p=i;
     }
  }

  return p;
}

/*************************************************************************
@
@-------------------------------------------------------------------------
@
@ static matrix_TYP **get_idem2(matrix_TYP **gen,int number,int *anz)
@
@
@-------------------------------------------------------------------------
@
**************************************************************************/
static matrix_TYP **get_idem2(matrix_TYP **gen,int number,int *anz)
{
   int i,
       j,
       k,
       found,
       dim_new1,
       dim_new2,
       vec[100];     /* I don't belivieve we will ever get a
                        centre of a centralizer which has bigger dimension */

   matrix_TYP *ONE,
              *e1,
              *e2,
              *tmp,
              *min,
              *zero,
             **gentrans,
             **gen_new1,
             **gen_new2,
             **erg;

   /* number == 0 is a technical statement */
   if (number == 0){
      anz[0] = 0;
      return NULL;
   }

   /* if the algbra has only got dimension one, we are finished */
   if (number == 1){
      min = min_pol(gen[0]);
      zero = zeros(min);
      erg = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP *));
      erg[0] = copy_mat(gen[0]);
      if (zero->cols == 1)
         erg[0]->kgv = zero->array.SZ[0][0];
      else
         erg[0]->kgv = zero->array.SZ[0][1];
      Check_mat(erg[0]);
      anz[0] = 1;
      free_mat(min);
      free_mat(zero);
      return erg;
   }

   /* calculate the right regular representation */
   /* assume the matrices in gen to be a Z-BASIS */
   gentrans = (matrix_TYP **) malloc(number * sizeof(matrix_TYP *));
   for (i=0;i<number;i++){
      gentrans[i] = init_mat(number,number,"");
      for (j=0;j<number;j++){
         tmp = mat_mul(gen[i],gen[j]);
         form_to_vec_modular(gentrans[i]->array.SZ[j],tmp,gen,number);
         free_mat(tmp);
      }
   }

   /* search for a matrix which has a reducible minimal polynomial */
   found = FALSE;
   for (i=0;i<number && !found;i++){
      min = min_pol(gentrans[i]);

      /* here is one part which has to be altered if one want's to
         be able to handle bigger fields, ie. do not only check whether
         min has got zero, but other factors */
      zero = zeros(min);
      if (zero->cols > 1 ||
         (zero->cols == 1 && min->cols > 3)){
         found = TRUE;
      }
      else{
         free_mat(min);
         free_mat(zero);
      }
   }
   i--;

   /* calculate the identity in this algebra (which sounds a bit crazy) */
   tmp = init_mat(number,number,"1");
   form_to_vec(vec,tmp,gentrans,number,&j);
   ONE = vec_to_form(vec,gen,number);
   ONE->kgv = j;
   free_mat(tmp);

   /* free gentrans, which isn't needed anymore */
   for (j=0;j<number;j++){
      free_mat(gentrans[j]);
   }
   free(gentrans);

   if (found){

      /* the definition of e1 is the next part which has to be altered
         formaly e1 = min/factor (gen[i]). if these two things are changed,
         we are able to handle bigger fields */

      /* split the algebra in two of it's components
         by multiplying whith the two factors of min */
      e1 = imat_add(gen[i],ONE,1,-zero->array.SZ[0][0]);
      div_by_root(min,zero->array.SZ[0][0]);

      gen_new1 = (matrix_TYP **) malloc(number * sizeof(matrix_TYP *));
      for (j=0;j<number;j++){
         gen_new1[j] = mat_mul(gen[j],e1);
      }
      dim_new1 = long_rein_formspace(gen_new1,number,0);

      /* free space before doing the recursion */
      free_mat(zero);
      free_mat(min);
      free_mat(e1);

      /* do the  first recursion */
      erg = get_idem2(gen_new1,dim_new1,&i);

      /* now ONE - \sum erg[i] will be an idempotent */
      e2 = copy_mat(ONE);
      for (j=0;j<i;j++){
         tmp = imat_add(e2,erg[j],1,-1);
         free_mat(e2);
         e2 = tmp;
      }

      gen_new2 = (matrix_TYP **) malloc(number * sizeof(matrix_TYP *));
      for (j=0;j<number;j++){
         gen_new2[j] = mat_mul(gen[j],e2);
         gen_new2[j]->kgv = 1;
         Check_mat(gen_new2[j]);
      }
      free_mat(e2);
      dim_new2 = long_rein_formspace(gen_new2,number,0);

      /* do the second recursion */
      gentrans = get_idem2(gen_new2,dim_new2,&j);

      /* copy  the idempotents to the right place */
      if (j>0) erg = (matrix_TYP **) realloc(erg,(i+j) * sizeof(matrix_TYP *));
      anz[0] = i+j;
      for (k=0;k<j;k++) erg[i+k] = gentrans[k];

      /* free the rest of the space */
      for (i=0;i<number;i++){
         free_mat(gen_new1[i]);
         free_mat(gen_new2[i]);
      }
      free(gen_new1);
      free(gen_new2);
      if (j>0) free(gentrans);

      /* remove duplicates in erg */
      for (i=1;i<anz[0];i++){
        if (pos(erg,i,erg[i]) != -1){
          free_mat(erg[i]);
          anz[0]--;
          erg[i] = erg[anz[0]];
        }
      }
   }
   else{
      erg = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP));
      anz[0] = 1;
      erg[0] = ONE;
   }

   if (found) free_mat(ONE);

   return erg;
}

void reduce_by_idem(matrix_TYP *id,matrix_TYP **cen,int no,int offset)
{
  int i,
      j,
      k;

  rational one,
           minus_one;

  matrix_TYP *lines,
             *tmp;

  one.n = 1;
  one.z =1;
  minus_one.n =1;
  minus_one.z =-1;

  for (i=offset;i<no;i++){
     put_mat(cen[i],NULL,"red",2);
     tmp = mat_mul(id,cen[i]);
     mat_addeq(cen[i],tmp,one,minus_one);
     Check_mat(cen[i]);
     cen[i]->kgv = 1;
     put_mat(cen[i],NULL,"red",2);
     free_mat(tmp);
  }

  lines = init_mat(no-offset,id->rows*id->cols,"");
  for (i=offset;i<no;i++)
     for (j=0;j<id->rows;j++)
        for (k=0;k<id->cols;k++)
           lines->array.SZ[i-offset][j*id->cols+k] = cen[i]->array.SZ[j][k];

  tmp = long_rein_mat(lines);

  free_mat(cen[offset]);
  cen[offset] = copy_mat(id);

  for (i=offset+1;i<no;i++)
     for (j=0;j<id->rows;j++)
        for (k=0;k<id->cols;k++)
           cen[i]->array.SZ[j][k] = tmp->array.SZ[i-1-offset][j*id->cols + k];

  for (i=0;i<no;i++) Check_mat(cen[i]);

  if (INFO_LEVEL & 4){
     put_mat(lines,NULL,"lines",2);
     put_mat(tmp,NULL,"tmp",2);
     for (i=0;i<no;i++){
        put_mat(cen[i],NULL,"cen in reduce_by_idem",2);
     }
  }

  free_mat(tmp);
  free_mat(lines);

  return;
} /* reduce_by_idem(...) */

/************************************************************************
@
@------------------------------------------------------------------------
@
@ matrix_TYP **idempotente(matrix_TYP **gen,int gen_no,matrix_TYP *form,
@                          int *anz,int *dimc,int *dimcc,int *options)
@
@------------------------------------------------------------------------
@
*************************************************************************/
matrix_TYP **idempotente(matrix_TYP **gen,int gen_no,matrix_TYP *form,
                         int *anz,int *dimc,int *dimcc,int *options)
{

  matrix_TYP **centralizer,   /* holds a Q-basis for the centralizer of
                                 QG (as enveloping algebra) */
             **ccentralizer,  /* holds a Q-basis for the centre of
                                 centralizer */
              *tmp,
              *tmp2,
              *test,
             **eigenspace,
              *min,           /* the minpol of a result to be */
              *forminv,       /* inverse of form */
              *action,        /* describes the action of form on
                                 ccentralizer via (form * Z * form^(-1))^Tr */
             **idem;          /* the result will be stucked in here */

  int i,
      j,
     *l,       /* holds a temporaly vector for vec_to_form */
      k,
      d,
      den,
      flag;

  /* transform all generators to kgv format */
  for (i=0;i<gen_no;i++){
     rat2kgv(gen[i]);
  }

  /* calculate the centralizer */
  centralizer = solve_endo(gen,gen,gen_no,dimc);

  /* calculate the centre of the centralizer */
  ccentralizer = calc_ccentralizer(centralizer,dimc[0],gen,gen_no,dimcc);

  /* output for debugging */
  if (INFO_LEVEL & 4){
     for (i=0;i<dimc[0];i++) put_mat(centralizer[i],NULL,"cen",2);
     for (i=0;i<dimcc[0];i++) put_mat(ccentralizer[i],NULL,"ccen",2);
  }

  /* calculate action, actualy d * action, because I don't like denominators.
     we calculate action as acting on rows, because this is very convient */
  forminv = long_mat_inv(form);
  rat2kgv(forminv);
  d = forminv->kgv;
  forminv->kgv=1;
  action = init_mat(dimcc[0],dimcc[0],"");

  for (i=0;i<dimcc[0];i++){
     tmp = mat_kon(form,ccentralizer[i],forminv);
     tmp2 = tr_pose(tmp);
     free_mat(tmp);
     tmp = tmp2;

     if (INFO_LEVEL & 4){
        for (j=0;j<dimcc[0];j++) put_mat(ccentralizer[j],NULL,"ccen",2);
        put_mat(tmp,NULL,"tmp",2);
        put_mat(form,NULL,"form",2);
        put_mat(forminv,NULL,"forminv",2);
     }

     /* express tmp as linear combination of the matrices in
        ccentralizer. This is possible in an integral way, because
        we have a Z-basis in ccentralizer */
     form_to_vec(action->array.SZ[i],tmp,ccentralizer,dimcc[0],&den);
     if (den != 1){
        fprintf(stderr,"error in idempotente, expression should be possible\n");
        fprintf(stderr,"in an integral way\n");
        exit(3);
     }
     free_mat(tmp);
  }

  /* calculate the 1 eigenspace of action, or strictly speeking the
     d eigenspace of d*action */
  for (i=0;i<action->rows;i++) action->array.SZ[i][i] -= d;
  tmp = tr_pose(action);
  free_mat(action);
  action = tmp;
  eigenspace = long_solve_mat(action,NULL);
  if (eigenspace[0] != NULL || eigenspace[1] == NULL){
     fprintf(stderr,"error in idempotente\n");
     exit(3);
  }
  tmp = eigenspace[1];
  free(eigenspace);
  eigenspace = (matrix_TYP **) malloc(tmp->cols*sizeof(matrix_TYP *));
  l = (int *) malloc(dimcc[0] * sizeof(int));
  for (i=0;i<tmp->cols;i++){
     for (j=0;j<tmp->rows;j++) l[j] = tmp->array.SZ[j][i];
     eigenspace[i] = vec_to_form(l,ccentralizer,dimcc[0]);
  }

  d = tmp->cols;

  if (INFO_LEVEL & 4){
     put_mat(action,NULL,"action",2);
  }

  /* clean up again */
  free(l);
  free_mat(tmp);
  free_mat(action);
  free_mat(forminv);

  if (INFO_LEVEL & 4){
    for (i=0;i<d;i++) put_mat(eigenspace[i],NULL,"eigenspace",2);
  }

  idem = get_idem2(eigenspace,d,anz);

  /* clean up the last space */
  for (i=0;i<d;i++) free_mat(eigenspace[i]);
  free(eigenspace);

  /* put the centralizer and it's centre on the same pointer as idem */
  idem = (matrix_TYP **) realloc(idem,(anz[0]+dimc[0]+dimcc[0])
                                     * sizeof(matrix_TYP));
  for (i=0;i<dimc[0];i++) idem[anz[0]+i] = centralizer[i];
  for (i=0;i<dimcc[0];i++) idem[anz[0]+dimc[0]+i] = ccentralizer[i];

  /* clear the pointer's needed */
  free(centralizer);
  free(ccentralizer);

  return idem;

}
