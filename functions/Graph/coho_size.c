#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <longtools.h>



/************************************************************************
@   Calculate the order of the cohomology group!
@
@   matrix_TYP *D:         2nd matrix returned by cohomolgy for G.
*************************************************************************/
MP_INT cohomology_size(matrix_TYP *D)
{
   int first, i;

   MP_INT coho_size;


   for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);

   mpz_init_set_si(&coho_size, 1);
   /* the order of the cohomology group in the product of those elementary
      divisors which are not equal 0 */
   for (i = first; i < D->cols && D->array.SZ[i][i] != 0; i++)
      mpz_mul_ui(&coho_size, &coho_size, (unsigned long) D->array.SZ[i][i]);

   return(coho_size);
}
