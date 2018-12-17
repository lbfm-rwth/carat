#include"typedef.h"
#include"matrix.h"
#include"symm.h"
#include"reduction.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: short_reduce.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *short_reduce(A, SV, Trf)
@ matrix_TYP *A, *SV, *Trf;
@
@ short_reduce make a reduction of the matrix A using the shortvectors
@ given in SV.
@ If an entry in SV is 1 or -1 the vector is used to make A better
@ afterwards the vectors in SV are transformed and the same is tried
@ again.
@ A and SV are not changed in this function.
@ the result is the reduced matrix and the transformation is
@ applied to the matrix Trf.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *short_reduce(A, SV, Trf)
matrix_TYP *A, *SV, *Trf;
{

   int step, i,j,k,l, dim;
   int min, min1, max, anz;
   int **B, **Titr ,**T, **S ;
   matrix_TYP *Aneu;
   int abbruch, lastchange, merk, *mp;
   int *v;

   dim = A->cols;
   min = SV->array.SZ[0][dim];
   anz = SV->rows;
   if( (v = (int *)malloc(dim *sizeof(int ))) == NULL){
     printf("malloc of 'v' in short_reduce failed\n");
     exit(2);
   }
   Aneu = copy_mat(A);
   B = Aneu->array.SZ;
   if(Trf != NULL)
     T = Trf->array.SZ;
   else
     T = NULL;
   S = SV->array.SZ;
   if( (Titr = (int **)malloc(dim *sizeof(int *))) == NULL){
     printf("malloc of 'Titr' in short_reduce failed\n");
     exit(2);
   }
   for(i=0;i<dim;i++)
   {
     if( (Titr[i] = (int *)malloc(dim *sizeof(int))) == NULL)
     {
       printf("malloc of 'Titr[%d]' in short_reduce failed\n", i);
       exit(2);
     }
     for(j=0;j<dim;j++)
       Titr[i][j] = 0;
     Titr[i][i] = 1;
   }
   min1 = B[0][0];
   max = B[0][0];
   for(i=1; i<dim;i++)
   {
     if(B[i][i] < min1)
        min1 = B[i][i];
     if(B[i][i] > max)
        max = B[i][i];
   }
   if(max == min)
   { for(i=0;i<dim;i++)
       free(Titr[i]);
       free(Titr);
       return(Aneu);
   }
   if(min1 == min)
   {
     abbruch = FALSE;
     for(step=0;step<dim && abbruch == FALSE;step++)
     {
       if(B[step][step] != min)
       {
          for(i=step+1;i<dim && B[i][i] != min ;i++);
          if(i == dim)
            abbruch = TRUE;
          else
          {
             if(T != NULL)
             {
               mp = T[step];
               T[step] = T[i];
               T[i] = mp;
             }
             mp = Titr[step];
             Titr[step] = Titr[i];
             Titr[i] = mp;
             mp = B[step];
             B[step] = B[i];
             B[i] = mp;
             for(j=0;j<dim;j++)
             {
                merk = B[j][step];
                B[j][step] = B[j][i];
                B[j][i] = merk;
             }
          }
       }
     }
   }
   for(step = 0;step<dim && B[step][step] == min; step++);
   lastchange = anz;
   abbruch = FALSE;
   while(step < dim && abbruch == FALSE)
   {
      for(i=0;i<anz && step < dim && i != lastchange ;i++)
      {
        /*----------------------------------------------------*\
        | Calculate S[i] * Titr^{tr}
        \*----------------------------------------------------*/
         for(k=step;k<dim;k++)
         {
           v[k] = 0;
           for(l=0;l<dim;l++)
             v[k] += S[i][l] * Titr[k][l];
         }
         for(k=step;k<dim && v[k] != 1 && v[k] != -1; k++);
         if(k<dim)
         {
           lastchange = i;
           for(j = 0;j<step;j++)
           {
             v[j] = 0;
             for(l=0;l<dim;l++)
               v[j] += S[i][l] * Titr[j][l];
           }
           if(v[k] == -1)
           {
             for(j=0;j<dim;j++)
               v[j] = -v[j];
           }
           if(k != step)
           {
              if(T != NULL)
              {
                mp = T[step];
                T[step] = T[k];
                T[k] = mp;
              }
              mp = Titr[step];
              Titr[step] = Titr[k];
              Titr[k] = mp;
              mp = B[step];
              B[step] = B[k];
              B[k] = mp;
              for(l=0;l<dim;l++)
              { merk = B[l][step]; B[l][step] = B[l][k]; B[l][k] = merk;}
              merk = v[step]; v[step] = v[k]; v[k] = merk;
           }
           for(j=0;j<step;j++)
           {
             if(v[j] != 0)
             {
                for(l=0;l<dim;l++)
                {
                  if(T != NULL)
                    T[step][l] += T[j][l] * v[j];
                  Titr[j][l] -= Titr[step][l] * v[j];
                  B[step][l] += B[j][l] * v[j];
                }
             }
           }
           for(j=step+1;j<dim;j++)
           {
             if(v[j] != 0)
             {
                for(l=0;l<dim;l++)
                {
                  if(T != NULL)
                    T[step][l] += T[j][l] * v[j];
                  Titr[j][l] -= Titr[step][l] * v[j];
/*
                  Titr[l][step] -= Titr[l][j] * v[j];
*/
                  B[step][l] += B[j][l] * v[j];
                }
             }
           }
           for(j=0;j<step;j++)
           {
             if(v[j] != 0)
             {
                for(l=0;l<dim;l++)
                  B[l][step] += B[l][j] * v[j];
             }
           }
           for(j=step+1;j<dim;j++)
           {
             if(v[j] != 0)
             {
                for(l=0;l<dim;l++)
                  B[l][step] += B[l][j] * v[j];
             }
           }
         step++;
         }
      }
      if(i == lastchange)
         abbruch = TRUE;
   }
   for(i=0;i<dim;i++)
     free(Titr[i]);
   free(Titr);
   free(v);
   return(Aneu);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pr_short_red(A, Trf)
@ matrix_TYP *A, *Trf;
@
@ The same as short_reduce but before and after using this function
@ a pair_redduction is used.
@ The shortest vectors are calculated by the function itsself.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *pr_short_red(A, Trf)
matrix_TYP *A, *Trf;
{
   matrix_TYP *Aneu, *Aneu1, *Aneu2, *SV;
   int i,j,k, min, max, min1, dim;
   int anz, len;
   int is_even;

   dim = A->cols;
   Aneu = copy_mat(A);
   is_even = TRUE;
   for(i=0;i<dim && is_even == TRUE; i++)
   {
      if( (A->array.SZ[i][i]%2) != 0)
        is_even = FALSE;
   }
   pr_red(Aneu->array.SZ, Trf->array.SZ, dim);
   min1 = Aneu->array.SZ[0][0];
   max = Aneu->array.SZ[0][0];
   for(i=1;i<dim;i++)
   {
     if(Aneu->array.SZ[i][i] < min1)
       min1 = A->array.SZ[i][i];  
     if(Aneu->array.SZ[i][i] > max)
       max = A->array.SZ[i][i];  
   }
   if(is_even == TRUE)
     len = min1 - 2;
   else
     len = min1 -1;

   /* changed 17/1/97 from:
   short_vectors(Aneu, len, 0, 1, 1, &anz);
   to: */
   SV = short_vectors(Aneu, len, 0, 1, 1, &anz);
   free_mat(SV);

   if(anz == 0 && max == min1)
      return(Aneu);
   SV = shortest(Aneu, &min);
   Aneu1 = short_reduce(Aneu, SV, Trf);
   free_mat(SV);
   free_mat(Aneu);
   if(Aneu1->array.SZ[dim-1][dim-1] == min)
      return(Aneu1);
   pr_red(Aneu1->array.SZ, Trf->array.SZ, dim);
   return(Aneu1);
}
