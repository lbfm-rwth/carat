#include "typedef.h"
#include "gmp.h"
#include "matrix.h"

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ static matrix_TYP *reverse_valuation(MP_INT *val,matrix_TYP *D)
@
@ The value of val will not be changed !
@----------------------------------------------------------------------------
@
*****************************************************************************/
matrix_TYP *reverse_valuation(MP_INT *val,matrix_TYP *D)
{
   int i,
       first,
       last;

   matrix_TYP *erg;

   MP_INT tmp,
          copy;

   mpz_init(&tmp);
   mpz_init_set(&copy,val);

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);

   erg = init_mat(last-first,1,"");

   for (i=first;i<last;i++){
      mpz_mod_ui(&tmp,&copy,(unsigned long) D->array.SZ[i][i]);
      erg->array.SZ[i-first][0] = mpz_get_si(&tmp);
      mpz_sub_ui(&copy,&copy,(unsigned long) erg->array.SZ[i-first][0]);
      mpz_div_ui(&copy,&copy,(unsigned long) D->array.SZ[i][i]);
   }

   mpz_clear(&tmp);
   mpz_clear(&copy);

   return erg;
}
