#include"typedef.h"
#include"matrix.h"
#include"symm.h"
#include"longtools.h"
#include"tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: mink_red.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


static void Gramtrans(G, T, U, Uinv, vec, step)
matrix_TYP *G, *T;
int **U, **Uinv;
int *vec;
int step;
{
  int i,j,k,n;
  int **u, **uinv;
  int **merk1, **merk2;

  n = G->cols;
  u = (int **)malloc(n *sizeof(int *));
  uinv = (int **)malloc(n *sizeof(int *));
  merk1 = (int **)malloc(n *sizeof(int *));
  merk2 = (int **)malloc(n *sizeof(int *));
  for(i=0; i<n;i++)
  {
   u[i] = (int *)malloc(n *sizeof(int));
   uinv[i] = (int *)malloc(n *sizeof(int));
   merk1[i] = (int *)malloc(n *sizeof(int));
   merk2[i] = (int *)malloc(n *sizeof(int));
  }
  for(i=0; i<step; i++)
  {
    for(j=0; j<n; j++)
    {
      u[i][j] = 0;
      uinv[i][j] = 0;
    }
    u[i][i] = 1;
    uinv[i][i] = 1;
  }
  for(i=step; i<n; i++)
  {
    for(j=0; j<step; j++)
    {
      u[i][j] = 0;
      u[j][i] = 0;
      uinv[i][j] = 0;
      uinv[j][i] = 0;
    }
  }
  for(i=step; i<n; i++)
  {
    for(j=step; j<n;j++)
    {
       u[i][j] = U[i-step][j-step];
       uinv[i][j] = Uinv[i-step][j-step];
    }
  }
  for(i=0; i<(n-step); i++)
    for(j=0; j<step; j++)
       u[i+step][j] = -(vec[j])* (U[i][0]);
  for(i=0; i<step; i++)
     uinv[step][i] = vec[i];
  for(i=0; i<n; i++)
  {
     for(j=0; j<n; j++)
     {
       merk1[i][j] = 0;
       for(k=0; k<n; k++)
          merk1[i][j] += T->array.SZ[i][k] * u[k][j];
     }
  }
  for(i=0; i<n; i++)
    for(j=0;j<n;j++)
      T->array.SZ[i][j] = merk1[i][j];
  for(i=0; i<n; i++)
  {
     for(j=0; j<n; j++)
     {
       merk1[i][j] = 0;
       for(k=0; k<n; k++)
          merk1[i][j] += uinv[i][k] * G->array.SZ[k][j];
     }
  }
  for(i=0; i<n; i++)
  {
     for(j=0; j<n; j++)
     {
        merk2[i][j] = 0;
        for(k=0; k<n; k++)
          merk2[i][j] += merk1[i][k] * uinv[j][k];
     }
  }
  for(i=0; i<n; i++)
    for(j=0;j<n;j++)
      G->array.SZ[i][j] = merk2[i][j];
  for(i=0; i<n; i++)
    { free(u[i]); free(uinv[i]); free(merk1[i]); free(merk2[i]);}
  free(u); free(uinv); free(merk1); free(merk2);
}


static void kurvectrans(SV, U, no, step)
matrix_TYP *SV;
int **U;
int no, step;
{
  int i,j, k;
  int n;
  int *vec;

  n = SV->cols;
  vec = (int *)malloc(n *sizeof(int));
  for(i=0;i<no;i++)
  {
    for(j=0; j<step; j++)
      vec[j] = SV->array.SZ[i][j];
    for(j=step;j<n; j++)
    {
      vec[j] = 0;
      for(k=step; k<n;k++)
        vec[j] += SV->array.SZ[i][k] * U[k-step][j-step];	
    }
    for(j=0; j<step; j++)
      vec[j] -= vec[step] * SV->array.SZ[no][j];
    for(j=0; j<n;j++)
      SV->array.SZ[i][j] = vec[j];
  }
  for(i=no+1;i<SV->rows;i++)
  {
    for(j=0; j<step; j++)
      vec[j] = SV->array.SZ[i][j];
    for(j=step;j<n; j++)
    {
      vec[j] = 0;
      for(k=step; k<n;k++)
        vec[j] += SV->array.SZ[i][k] * U[k-step][j-step];	
    }
    for(j=0; j<step; j++)
      vec[j] -= vec[step] * SV->array.SZ[no][j];
    for(j=0; j<n;j++)
      SV->array.SZ[i][j] = vec[j];
  }
  for(i=0; i<n; i++)
    SV->array.SZ[no][i] = 0;
  SV->array.SZ[no][step] = 1;
  free(vec);
}

static int **trafoinverse(vec, step, dim)
int *vec;
int step, dim;
{
  int i,j, s, k, v, w, udet, merk;
  int **U;
  int n;

  n = dim - step;
  U = (int **)malloc(n *sizeof(int *));
  for(i=0; i<n; i++)
   U[i] = (int *)malloc(n *sizeof(int));

  for(i=0; i<n;i++)
   U[0][i] = vec[i+step];
  for(k=0;k<n && U[0][k] == 0; k++);
  for(i=1; i<=k;i++)
   for(j=0;j<n;j++)
     U[i][j] = 0;
  for(i=1; i<=k;i++)
     U[i][i-1] = 1;
  udet = U[0][k];
  for(s=k+1; s<n;s++)
  {
    if(U[0][s] == 0)
    {
       for(i=0; i<s; i++)
        U[s][i] = 0;
       for(i=s+1; i<n; i++)
        U[s][i] = 0;
       U[s][s] = 1;
    }
    else
    {
       gcd_darstell(udet, U[0][s], &v, &w, &merk);
       U[s][s] = v;
       for(i=0; i<s; i++)
         U[s][i] = -(U[0][i]/udet) * w;
       udet = merk;
       for(i=s+1; i<n; i++)
        U[s][i] = 0;
    }
  }
  return(U);
}

static int wechsel(vec, step, udim)
int *vec;
int step, udim;
{
  int i, j;
  int a,b,c;
  i=step;
  while(i<udim && vec[i] == 0)
     i++;
  if(i== udim)
    return(FALSE);
  a = vec[i];
  for(j=i+1; j<udim && a != 1 && a!= -1; j++)
  {
    if(vec[j] != 0)
    {
      b = vec[j];
      while ((c=a%b)!=0){a=b; b=c;}
      a = b;
    }
  }
  if(a == 1 || a == -1)
    return(TRUE);
  return(FALSE);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mink_red(G, Trf)
@ matrix_TYP *G, *Trf;
@
@ calculates a matrices A and Trf such that A = Trf * G * Trf^{tr}
@ is Minkowski_reduced
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *mink_red(G, Trf)
matrix_TYP *G, *Trf;
{
  int i, j, k, l, step, anz;
  int n, bound, merk;
  int **U, **Uinv;
  int ergmax;
  matrix_TYP *erg;
  matrix_TYP *kurvecs;
  matrix_TYP *UU, *UUi;
  matrix_TYP *T, *Ti;

  extern int **trafoinverse();
  extern int wechsel();
  extern void Gramtrans();
  extern void kurvectrans();
  extern int matrixinverses();

  n = G->cols;

  /* ------------------------------------------------------------------- *\
  |     setze T = I_n  und erg = G                                        |
  \* ------------------------------------------------------------------- */
  T = init_mat(n,n,"1");
  erg = init_mat(G->rows, G->cols, "");
  for(i=0;i<G->rows;i++)
    for(j=0;j<G->cols;j++)
      erg->array.SZ[i][j] = G->array.SZ[i][j];
  /* ------------------------------------------------------------------- *\
  |     ordne erg tolal                                                     |
  \* ------------------------------------------------------------------- */
  for(i=0; i<n; i++)
  {
    k=i;
    for(j=i+1; j<n; j++)
    {
      if(erg->array.SZ[j][j] < erg->array.SZ[k][k])
        k=j;
    }
    if(k!=i)
    {
       for(j=0; j<n; j++)
       {
          l=erg->array.SZ[i][j];
          erg->array.SZ[i][j] = erg->array.SZ[k][j];
          erg->array.SZ[k][j] = l;
       }
       for(j=0; j<n; j++)
       {
          l=erg->array.SZ[j][i];
          erg->array.SZ[j][i] = erg->array.SZ[j][k];
          erg->array.SZ[j][k] = l;
          l=T->array.SZ[j][i];
          T->array.SZ[j][i] = T->array.SZ[j][k];
          T->array.SZ[j][k] = l;
       }
    }
  }

  /* ------------------------------------------------------------------- *\
  |     Allocieren und Initialisierung des Speicherplatzes                |
  \* ------------------------------------------------------------------- */
  erg->flags.Symmetric = TRUE;
  if(n == 0 || n == 1)
     return(erg);

  if((UUi = (matrix_TYP *)malloc(sizeof(matrix_TYP))) == NULL){
    printf("malloc of 'UUi' in 'mink_red' failed\n");
    exit(2);
  }
  UUi->kgv = 1;
  ergmax = erg->array.SZ[n-1][n-1];
  kurvecs = short_vectors(erg, ergmax, 0, 0,0,&anz);

  for(step = 0; step < n;step++)
  {
    bound = erg->array.SZ[step][step];
    for(i=0; i<kurvecs->rows; i++)
    {
     if(kurvecs->array.SZ[i][n] < bound && wechsel(kurvecs->array.SZ[i],step,n)== 1)
     {
        Uinv = trafoinverse(kurvecs->array.SZ[i], step, n);
        UUi->cols = UUi->rows = n-step;
        UUi->array.SZ = Uinv;
        UU = long_mat_inv(UUi);
        U = UU->array.SZ;
        Gramtrans(erg, T, U, Uinv, kurvecs->array.SZ[i], step);
        kurvectrans(kurvecs, U, i, step);
        bound = erg->array.SZ[step][step];
        for(j=0; j<(n-step); j++)
          free(Uinv[j]);
        free_mat(UU); free(Uinv);
     }
    }
  }

  /* printf("step = %d\n", step);
  put_mat(erg, NULL, "erg in mink_red", 2); */

  for(i=1; i<n;i++)
  {
    if(erg->array.SZ[i][i-1] < 0)
    {
      for(j=0; j<n; j++)
      {
         erg->array.SZ[i][j] = - erg->array.SZ[i][j];
         erg->array.SZ[j][i] = - erg->array.SZ[j][i];
         T->array.SZ[j][i] = - T->array.SZ[j][i];
      }
    }
  }
  erg->kgv =  G->kgv;
  if(erg->kgv != 1)
    erg->flags.Integral = FALSE;

  /* changed 16/1/97 tilman from
  free(kurvecs);
  to: */
  free_mat(kurvecs);

  free(UUi);
  Ti = long_mat_inv(T);
  free_mat(T);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
       Trf->array.SZ[i][j] = Ti->array.SZ[i][j];
  free_mat(Ti);

  /* put_mat(erg, NULL, "erg in mink_red", 2); */

  return(erg);
}
