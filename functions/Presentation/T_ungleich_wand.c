
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: T_ungleich_wand.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

/**************************************************************************\
@--------------------------------------------------------------------------
@ wall_TYP **T_ungleich_wand(vec,mat,grn,length)
@
@ matrix_TYP *mat:   the rows of mat span a lattice T
@ matrix_TYP *vec:   affine vector x, the starting point.
@ matrix_TYP *grn:   grammatrix 
@ int *length:       will denote the number of inequalities
@
@ calculates the inenqualities describing the faces of the Dirichletdomain
@ D(x, Tx)
@ 2tAx + tAt -2*tAy > 0  (t a row of mat, x =vec, A = grn)
@
@--------------------------------------------------------------------------
\**************************************************************************/
/* Dieses Programm berechnet die Ungleichungen des Dirichletbereichs
   D(x,Tx) eines Punktes x unter dem Translationsnormalteiler einer
   Raumgruppe R.
   Eingabe: Vektor x, relevante Translationen (als Matrix mat)
   und Grammatrix A des zugehoerigen Skalarproduktes (A=grn=form).
   ((Beachte: Punktgruppe ist aus der GL(n,Z) wenn das zugrundeliegende
     Skalarprod. gerade grn ist.))
   Ungleichung:   2tAx + tAt -2*tAy > 0  (t relevante Translation, x Startwert)
*/


#include"typedef.h"
#include"matrix.h"
#include"polyeder.h"
#include"getput.h"
#include"presentation.h"

wall_TYP **T_ungleich_wand(vec,mat,grn,length)
matrix_TYP *mat, *vec, *grn;
int *length;
{
  int i,j,anz;
  matrix_TYP *hmat,*out, *erg1, *erg,*tmp1, *tmp2, *tmp, *konst, *konst1;
  wall_TYP **wall;
  
  extern GGT();
  
  anz = 2*mat->rows; 
  *length = anz;
  wall = (wall_TYP **)malloc(anz *sizeof(wall_TYP *));
  
  for(i=0;i<anz;i++)
    wall[i] = init_wall(mat->cols+1);
  
  tmp1 = copy_mat(mat);
  sscal_mul(tmp1, -1);
  Check_mat(tmp1);
  real_mat(tmp1, tmp1->rows+(mat->rows), tmp1->cols);
  
  for(i = mat->rows; i < tmp1->rows; i++)
    for(j = 0; j < mat-> cols; j++)
      tmp1->array.SZ[i][j] = mat->array.SZ[i-mat->rows][j];
  
  tmp2 = mat_mul(grn,vec); Check_mat(tmp2);
  out = mat_mul(tmp1,tmp2);
  sscal_mul(out, 2);
  Check_mat(out);			/* out = 2 tAx */
  tmp = tr_pose(tmp1);
  konst1 = mat_mul(grn,tmp); Check_mat(konst1);
  konst = mat_mul(tmp1,konst1); Check_mat(konst);	/* konst = tAt */
  
  for(i = 0; i < out->rows; i++){     
    out->array.SZ[i][0]= (konst->kgv)*out->array.SZ[i][0]; 
    out->array.SZ[i][0] += (out->kgv)*konst->array.SZ[i][i];
  }
  out->kgv = (out->kgv)*(konst->kgv);
  Check_mat(out);			/* out = 2 tAx + tAt */
  free_mat(konst); konst = NULL;
  erg1 = copy_mat(tmp1);
  erg = mat_mul(erg1,grn);
  sscal_mul(erg, -2);        Check_mat(erg);	/* erg = (-2tA) */
  sscal_mul(erg, out->kgv);  Check_mat(erg);
  real_mat(erg, erg->rows, erg->cols+1);
  
  for(i = 0; i < erg->rows; i++)
    erg->array.SZ[i][erg->cols-1] = out->array.SZ[i][0];
  erg->kgv = (erg->kgv)*(out->kgv);
  Check_mat(erg);			/* erg = (-2tA, 2 tAx + tAt) */
  
  for(i=0;i<anz;i++)
    {
      wall[i]->dim = mat->cols+1;
      hmat = init_mat(mat->cols+1,mat->cols+1,"1");
      wall[i]->mat = hmat; 
      
      for(j= 0;j<erg->cols;j++)
	wall[i]->gl[j] = erg->array.SZ[i][j];
      normal_wall(wall[i]); 
      for(j=0;j<mat->cols;j++)
	wall[i]->mat->array.SZ[j][mat->cols] = tmp1->array.SZ[i][j];
      Check_mat(wall[i]->mat);
    }
  
  free_mat(out); out = NULL;
  free_mat(erg1); erg1 = NULL;
  free_mat(tmp); tmp = NULL;
  free_mat(tmp1); tmp1 = NULL; /* anne 3.9. */
  free_mat(tmp2); tmp2 = NULL;
  free_mat(konst1); konst1 = NULL;
  free_mat(erg); erg = NULL;
  
  return(wall);
}
