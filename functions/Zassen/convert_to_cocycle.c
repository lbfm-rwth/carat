#include "typedef.h"
#include "matrix.h"
#include "gmp.h"
#include "zass.h"

/**************************************************************************
@
@ -------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP *convert_to_cocycle(matrix_TYP *x,matrix_TYP *cocycle,
                               matrix_TYP *D)
{

   int i,
       j,
       first,
       last,
       factor;

   matrix_TYP *erg;

   erg = init_mat(cocycle->rows,1,"");

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);

   for (i=first;i<last;i++){
      factor = x->array.SZ[i-first][0] *
               D->array.SZ[last-1][last-1] / D->array.SZ[i][i];
      for (j=0;j<erg->rows;j++){
         erg->array.SZ[j][0] += (cocycle->array.SZ[j][i-first] * factor);
      }
   }

   erg->kgv = D->array.SZ[last-1][last-1];
   erg->flags.Integral = FALSE;

   Check_mat(erg);

   return erg;
}

