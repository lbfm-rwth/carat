#include "typedef.h"
#include "matrix.h"
#include "getput.h"
#include "datei.h"
#include "longtools.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: gittstab.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

static int 
is_in (matrix_TYP *M, bravais_TYP *B)
{
  int i,j,k;
  matrix_TYP *mtr;
  matrix_TYP *F, *A;
  int erg;
  erg = TRUE;

  mtr = init_mat(M->cols, M->rows, "");
  for(i=0; i<M->cols; i++)
    for(j=0; j<M->rows; j++)
      mtr->array.SZ[i][j] = M->array.SZ[j][i];

  for(i=0; i<B->form_no && erg == TRUE; i++)
  {
    A = mat_mul(mtr, B->form[i]);
    F = mat_mul(A, M);
    free_mat(A);
    for(j=0; j<F->rows && erg == TRUE; j++)
      for(k=0; k<F->cols && erg == TRUE; k++)
        if(F->array.SZ[j][k] != B->form[i]->array.SZ[j][k])
           erg = FALSE;
    free_mat(F);
  }
  free_mat(mtr);
  return(erg);
}



/*--------------------------------------------------------------------*\
|  test if a matrix Y is already contained in a list of matices        |
|  with length 'anz'.                                                  |
|  If Y is not contained in the list, the result is -1, else           |
|  the result is the position where Y was found.                       |
\*--------------------------------------------------------------------*/
static int 
such (matrix_TYP *Y, matrix_TYP **linv, int anz)
{
   int schonda, i;
   matrix_TYP *G;

   schonda = -1;
   for(i=0; i<anz && schonda == -1; i++)
   {
      G = mat_mul(linv[i], Y);
      Check_mat(G);
      if(G->kgv == 1)
        schonda = i;
      free_mat(G);
   }
   return schonda;
}

/*--------------------------------------------------------------------*\
| if A is not the identity-matrix and A is not contained in the        |
| geneator-matrices of stab the result is TRUE, else FALSE             |
\*--------------------------------------------------------------------*/
static int 
mat_such (bravais_TYP *stab, int **A)
{
   int i,j,k, neu;
         i=0; neu = 0;
         while(i<stab->dim && neu == 0)
         {
              j=0;
              while(j<stab->dim && neu == 0)
              {
                 if(i==j && A[i][j] != 1)
                    neu = 1;
                 if(i != j && A[i][j] != 0)
                    neu = 1;
                 j++;
              }
              i++;
         }
         if(neu == 0)
            return(neu);

   if(stab->gen_no == 0)
      return(1);
   neu = 1;
   i=0;
   do
   {
     j=0; k = stab->dim;
     while(j<stab->dim && k == stab->dim)
     {
        k=0;
        while(k<stab->dim && A[j][k] == stab->gen[i]->array.SZ[j][k])
          k++;
        j++;
     }
     if(j == stab->dim && k == stab->dim)
        neu =0;
     i++;
   }while (neu == 1 && i<stab->gen_no);
   return(neu);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *Z_class(B, zen)
@ bravais_TYP *B;
@ matrix_TYP *zen;
@
@---------------------------------------------------------------------------
@
@ calculates the intersection of zen^(-1)*B*zen with GL_n(Z)           |
\**************************************************************************/
bravais_TYP *
Z_class (bravais_TYP *B, matrix_TYP *zen)
{
   int i,
       j;
   bravais_TYP *S;
   bravais_TYP *C;
   bravais_TYP *N;
   bravais_TYP *CS = 0;
   bravais_TYP *NS;
   matrix_TYP *A1, *A2;
   matrix_TYP *zinv,
              *ztr;
   int anz;
   int *noetig = NULL;

   S = gittstabneu(B, zen);

   /* deleted this output 25/4/97 tilman
   put_bravais(S, NULL, "Stabilisator noch nicht konjugiert"); */

   zinv = mat_inv(zen);
   for(i=0; i<S->gen_no; i++)
   {
      A1 = mat_mul(zinv, S->gen[i]);
      A2 = mat_mul(A1, zen);
      free_mat(A1); free_mat(S->gen[i]);
      Check_mat(A2);
      S->gen[i] = A2;
   }
   C = (bravais_TYP *)calloc(1, sizeof(bravais_TYP) );
  C->gen_no = B->cen_no;
  C->dim = B->dim;
    C->gen = (matrix_TYP **) malloc(C->gen_no *sizeof(matrix_TYP *));
  for(i=0; i<C->gen_no; i++)
     C->gen[i] = B->cen[i];
  if(C->gen_no > 0)
  CS = gittstabneu(C, zen);
   S->form = (matrix_TYP **) malloc(B->form_no *sizeof(matrix_TYP *));
   S->form_no = B->form_no;
   ztr = init_mat(zen->cols, zen->rows, "");
   for(i=0; i<zen->rows; i++)
     for(j=0; j<zen->cols; j++)
       ztr->array.SZ[j][i] = zen->array.SZ[i][j];
   for(i=0; i<B->form_no; i++)
   {
      A1 = mat_mul(ztr, B->form[i]);
      S->form[i] = mat_mul(A1, zen);
      Check_mat(S->form[i]);
      S->form[i]->kgv = 1;
      S->form[i]->flags.Integral = TRUE;
      free_mat(A1);
   }

   /* inserted 05/05/97 tilman: do an rein on the formspace of S
      to give a Z-basis */
   long_rein_formspace(S->form,S->form_no,1);

   if(C->gen_no > 0)
   {
   S->cen_no = CS->gen_no;
   S->cen = (matrix_TYP **) malloc(S->cen_no *sizeof(matrix_TYP *));
   for(i=0; i<CS->gen_no; i++)
   {
      A1 = mat_mul(zinv, CS->gen[i]);
      A2 = mat_mul(A1, zen);
      free_mat(A1);
      Check_mat(A2);
      S->cen[i] = A2;
   }
   }
 if(B->normal_no == 1)
   Check_mat(B->normal[0]);
 if(B->normal_no > 1 || (B->normal_no == 1 && B->normal[0]->flags.Scalar == FALSE))
 {
  N = (bravais_TYP *)calloc(1, sizeof(bravais_TYP));
  N->gen_no = B->normal_no + B->gen_no + B->cen_no;
  N->dim = B->dim;
    N->gen = (matrix_TYP **) malloc(N->gen_no *sizeof(matrix_TYP *));
  for(i=0; i<B->gen_no ; i++)
    N->gen[i] = B->gen[i];
  for(i=0; i<B->normal_no; i++)
    N->gen[i+B->gen_no] = B->normal[i];
  for(i=0; i<B->cen_no; i++)
    N->gen[i+B->gen_no+B->normal_no] = B->cen[i];
  NS = gittstabneu(N, zen);
  noetig = (int *) malloc(NS->gen_no *sizeof(int));
  for(i=0; i<NS->gen_no; i++)
  {
     noetig[i] = TRUE;
     if(is_in(NS->gen[i], B) == TRUE)
        noetig[i] = FALSE;
     if(noetig[i] == TRUE)
     {
        A1 = mat_inv(NS->gen[i]);
        for(j=0; j<i && noetig[i] == TRUE; j++)
        {
           A2 = mat_mul(A1, NS->gen[j]);
           if(is_in(A2, B) == TRUE)
               noetig[i] = FALSE;
           free_mat(A2);
        }
        free_mat(A1);
     }
  }
  anz = 0;
  for(i=0; i<NS->gen_no; i++)
    if(noetig[i] == TRUE)
       anz++;
   if(anz == 0)
     anz = 1;
   S->normal_no = anz;
   S->normal = (matrix_TYP **) malloc(S->normal_no *sizeof(matrix_TYP *));

   anz = 0;
   for(i=0; i<NS->gen_no; i++)
   {
	     if(noetig[i] == TRUE)
	     {
	      A1 = mat_mul(zinv, NS->gen[i]);
	      A2 = mat_mul(A1, zen);
	      free_mat(A1);
	      Check_mat(A2);
	      S->normal[anz] = A2;
	      anz++;
	     }
	   }
	   free_bravais(NS);
	   free(N->gen);
	   free(N);
	   if(anz == 0)
	     S->normal[0] = einheitsmatrix(B->dim);
	  }
	  if(B->normal_no == 1 && B->normal[0]->flags.Scalar == TRUE)
	  {
	     S->normal_no = 1;
	     S->normal = (matrix_TYP **) malloc(S->normal_no *sizeof(matrix_TYP *));
	     S->normal[0] = einheitsmatrix(B->dim);
	  }
	  if(C->gen_no > 0) {
	   free_bravais(CS);
	  }
	   free(C->gen);
	   free(C);
	   free_mat(zinv); free_mat(ztr);

	 /* inserted 7/6/97 tilman */
	 if (noetig != NULL) free(noetig);

 return(S);
}





/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *gittstab(grp, X)
@ bravais_TYP *grp;
@ matrix_TYP *X;
@
@   gittstab calculates the stabilizer of a lattice X under the        |
@   operation of the group 'grp' acting via left-multiplication        |
@---------------------------------------------------------------------------
@
\**************************************************************************/
bravais_TYP *
gittstab (bravais_TYP *grp, matrix_TYP *X)
{
   int i, j;
   matrix_TYP **list;
   matrix_TYP **list_inv;
   matrix_TYP **history;
   matrix_TYP *Y, *Yinv;
   matrix_TYP *Z;
   matrix_TYP *G, *Hinv;
   bravais_TYP *stab;
   int *index;
   int num, erz, schonda;
   int neu, anz = 0;

   if(grp->dim != X->rows)
   {
       printf("Matrixgruppe und Vektorraum haben unterschiedliche Dimension\n");
           exit(3);
   }
   if(X->cols != X->rows)
   {
       printf("Teilgitter hat nicht vollen Rang\n");
           exit(3);
   }
   list = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
   list_inv = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
   history = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
   stab = (bravais_TYP *) calloc(1, sizeof(bravais_TYP) );
   stab->dim = grp->dim;
   stab->gen_no = 0;
   stab->form_no = 0;
   stab->zentr_no = 0;
   stab->normal_no = 0;
   stab->order = 0;
   Z = init_mat(X->rows, X->rows, "");
   for(i=0; i<X->rows; i++)
     Z->array.SZ[i][i] = 1;
   
   Y = init_mat(X->rows, X->rows, "");
   for(i=0; i<X->rows; i++)
      for(j=0; j<X->cols; j++)
         Y->array.SZ[i][j] = X->array.SZ[i][j];
   Yinv = mat_inv(Y);
   list[anz] = Y;  /* init the orbit of X under 'grp' */
   list_inv[anz] = Yinv;
   history[anz] = Z; /* history[n] * list[0] = list[n] */
   anz++;
   num = 0;
   erz = 0;
/*--------------------------------------------------------------------*\
|  algorithmen to calculate the orbit of X                             |
\*--------------------------------------------------------------------*/
   do
   {
        Y = mat_mul(grp->gen[erz], list[num]);
        Z = mat_mul(grp->gen[erz], history[num]);
        schonda = such(Y, list_inv, anz);
        /*------------------------------------------------------------*\
        |  if new lattice in the orbit                                 |
        \*------------------------------------------------------------*/
        if (schonda == -1)
        {
          list = (matrix_TYP **) realloc(list, (anz+1) *sizeof(matrix_TYP *));
          list_inv = (matrix_TYP **) realloc(list_inv, (anz+1) *sizeof(matrix_TYP *));
          history = (matrix_TYP **) realloc(history, (anz+1) *sizeof(matrix_TYP *));
          Yinv = mat_inv(Y);
          list[anz] = Y;
          list_inv[anz] = Yinv;
          history[anz] = Z;
          anz++;
        }
        /*------------------------------------------------------------*\
        |  if lattice already found, calculate stabilizer-element      |
        \*------------------------------------------------------------*/
        else
        {
           Hinv = mat_inv(history[schonda]);
           G = mat_mul(Hinv, Z);
           neu = mat_such(stab, G->array.SZ);
           if(neu == 1)
           {
              if(stab->gen_no != 0)
                stab->gen = (matrix_TYP **) realloc(stab->gen, (stab->gen_no+1) *sizeof(matrix_TYP *));
              if(stab->gen_no == 0)
                stab->gen = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
              stab->gen[stab->gen_no] = G;
              stab->gen_no++;
           }
           else
            free_mat(G);
           free_mat(Hinv); free_mat(Z); free_mat(Y);
        }
        if (num == anz-1 && erz == grp->gen_no-1)
	        break;
        if (erz < grp->gen_no-1)
        ++erz;
        else
        {erz = 0;
        ++num;}
   }while (1);
/**************
   for(i=0; i<anz; i++)
     put_mat(list_inv[i], NULL, "list_inv" , 2);
*****************/
/*--------------------------------------------------------------------*\
| calculation of the order of the stabilizer                           |
\*--------------------------------------------------------------------*/
   index = factorize(anz);
   for(i=1; i<100; i++)
     stab->divisors[i] = grp->divisors[i] - index[i];
   if(grp->divisors[0] != 0 || index[0] != 0)
   {
     stab->divisors[0] = 1;
     stab->order = 0;
   }
   else
   {
     stab->order = 1;
     for(i=2; i<100; i++)
       for(j=0; j<stab->divisors[i]; j++)
         stab->order *= i;
   }
   for(i=0; i<anz; i++)
   {
      free_mat(list[i]);
      free_mat(list_inv[i]);
      free_mat(history[i]);
   }
   free(list);
   free(list_inv);
   free(history);
   free(index);
   return(stab);
}
