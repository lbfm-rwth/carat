/* berechne Grammatrix, hier bilde M*A*M^tr mit A=invariantes Skalarprodukt 
   (obiges gilt fuer Zeilenkonvention, bei Spalten bilde M^tr*A*M)          */

#include"typedef.h"
#include"matrix.h"

/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *neugram(mat, A)
@ matrix_TYP *mat, *A;
@
@ calculates M*A*M^tr 
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *neugram(mat, A)
matrix_TYP *mat, *A;

{
   matrix_TYP *erg, *tmp, *tmp1;

	extern matrix_TYP *mat_mul();
	extern matrix_TYP *tr_pose();
        tmp = tr_pose(mat);
	tmp1 = mat_mul(mat, A);
	erg = mat_mul(tmp1, tmp);
        free_mat(tmp); tmp = NULL;
        free_mat(tmp1); tmp1 = NULL;
	return(erg);
}
