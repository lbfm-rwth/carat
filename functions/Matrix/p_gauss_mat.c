#include "typedef.h"
#include "tools.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: p_gauss_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  p_gauss*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ int p_gauss (L_mat)
@ matrix_TYP *L_mat ;
@
@ applies a gauss algorithm to the rows of L_mat modulo L_mat->prime.
@ The entries of the matrix are changed.
@ The return is the rank of L_mat modulo L_mat->prime.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
p_gauss (matrix_TYP *L_mat)
{  
int **M, *v, f, help ;
int *M_j, *M_act, numax;
int *ss, *nuin;
int i, j, k, kk, cL, cM, rL, rM, rg, act;

  cM = cL = L_mat->cols;
  rM = rL = L_mat->rows;
  
  M   = L_mat->array.SZ;
  nuin= (int  *)malloc((1+cM)*sizeof(int));
  ss  = (int  *)malloc(cM*sizeof(int));
  
  rg = 0;
    /*------------------------------------------------------------*\
    | Gauss Algorithm                       |
    \*------------------------------------------------------------*/
  for ( i = 0; i < cM; i++ ) {
#if 0
    /*{{{  auskommentiert, printout*/
    for(pri = 0; pri < rM; pri++) {
      for(prj = 0; prj < cM; prj++) {
        if(M[pri][prj] != 0) {
          printf("%d ",M[pri][prj]);
        } else {
          printf("  ");
        }
      }
      printf("\n");
    }
    printf(" i = %d \n", i);
    /* ende der Einfuegung */
    /*}}}  */
#endif
    ss[i] = -1;
    /*---------------------------------------------------------*\
    | Find the row with non-zero entry in i-th column           |
    \*---------------------------------------------------------*/
    for ( j = rg; (j < rM) && (M[j][i] == 0); j++);
    if ( j == rM ) {
      continue;    
    }
    act = j;
    /*---------------------------------------------------------*\
    | Normalize act-th row and mark the non-zero entries        |
    \*---------------------------------------------------------*/
    nuin[0] = 1;
    f = P(1,-M[j][i]);
    for ( k = i ; k < cM; k++ ) {
      if((help = M[act][k])) {
        M[act][k] = P(help,f);
        nuin[nuin[0]++] = k;
      }                          
    }
    /*---------------------------------------------------------*\
    | Swap act-th row and rg-th row                  |
    \*---------------------------------------------------------*/
    v      = M[rg];
    M[rg]  = M[act];
    M[act] = v;
    ss[rg] = i;
    act    = rg ++;
    /*---------------------------------------------------------*\
    | Clear i-th column  downwards                |
    \*---------------------------------------------------------*/
    M_act = M[act];
    for ( j = act+1; j < rM; j++) {
      if ( (f = S(0,-M[j][i])) != 0 ) {
        M_j = M[j];
        M_j[i] = 0;
        numax = nuin[0];
        if( f == 1) {
          for(k=2; k < numax; ++k ) {
	    kk = nuin[k];
            M_j[kk] = S(M_j[kk],M_act[kk]);                
          }         
        } else {
          for(k=2; k < numax; ++k ) {
	    kk = nuin[k];
            M_j[kk] = S(M_j[kk],P(M_act[kk],f));           
          }
        }
      }
    }
  }
    /*========================================================*\
    || clear it upwards                                       ||
    \*========================================================*/
  for (i = rg-1; i > 0; i--) {
    nuin[0] = 2;
    nuin[1] = ss[i];
    M_act = M[i];
    for(j = ss[i]+1; j < cM; j++) {
      if(M_act[j]) {
        nuin[nuin[0]++] = j;
      }
    }
    if(nuin[0] == 2) {
      j = ss[i];
      for (k = i-1; k > 0; k--) {
        M[k][j] = 0;            
      }
    } else {
      for ( j = i-1; j >= 0; j--) {
        if ( (f = S(0,-M[j][ss[i]])) != 0 ) {
          M_j = M[j];
          M_j[ss[i]] = 0;
          numax = nuin[0];
          if( f == 1) {
            for(k=2; k < numax; ++k ) {
	      kk = nuin[k];
              M_j[kk] = S(M_j[kk],M_act[kk]);                 
            }         
          } else {
            for(k=2; k < numax; ++k ) {
	      kk = nuin[k];
              M_j[kk] = S(M_j[kk],P(M_act[kk],f));           
            }
          }
        }
      }
    }
  }
  free (ss);
  free (nuin);
  return (rg);
}
/*}}}  */
