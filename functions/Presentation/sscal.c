#include"typedef.h"
#include"matrix.h"

	/*============================================================*\
	|| Multiply matrix with short-integer scalar                  ||
	\*============================================================*/

/**************************************************************************\
@--------------------------------------------------------------------------
@ void sscal_mul(mat, v)
@ matrix_TYP *mat;
@ int v;
@
@ Multiply matrix with short-integer scalar  
@
@--------------------------------------------------------------------------
\**************************************************************************/
void sscal_mul(mat, v)
matrix_TYP *mat;
int v;
{
int i,j;
int **SZ;

/* handle the trivial cases */
if (v == 1) return;
if (mat->kgv == 0) return;
if (v == 0) {
	for(i = 0; i < mat->rows; i++) 
	for(j = 0; j < mat->cols; j++) 
		mat->array.SZ[i][j] = 0;
	return;
	}

SZ = mat->array.SZ;
if(mat->flags.Diagonal) {
	for(i = 0; i < mat->rows; i++) 
		SZ[i][i] *= v;
	}
else {
	for(i = 0; i < mat->rows; i++) 
	for(j = 0; j < mat->cols; j++) 
		SZ[i][j] *= v;
	}

}

