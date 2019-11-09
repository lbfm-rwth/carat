#include "typedef.h"
#include "matrix.h"
#include "tools.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: elt_div_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  local typedefs*/
/*
 *  local typedefs, only used in this file.
 */
typedef struct 
{ 
  int  f1, f2, g1, g2;   
} elt_div_pair;
typedef struct 
{ 
  int low, hi;
} elt_div_prod;

/*}}}  */
/*{{{  mindivmod, static*/
	/*============================================================*\
	||                                                            ||
	||  Loest azahl = prod.hi * bzahl + prod.low mit (meistens)   ||
	||  prod.low <= bzahl/2                                       ||
	||                                                            ||
	\*============================================================*/

elt_div_prod 
ZZ_mindivmod (int azahl, int bzahl)
{
elt_div_prod c;
int q_signum, r_signum;

	c.hi = c.low = 0;
	r_signum = azahl ? ( azahl > 0 ? 1 : -1 ) : 0;
	q_signum = bzahl ? ( bzahl > 0 ? 1 : -1 ) : 0;
	q_signum = r_signum == q_signum ? 1 : -1;

	azahl= abs(azahl);
        bzahl= abs(bzahl);
        c.hi = azahl / bzahl;
	c.low= azahl % bzahl;
	if ( (c.low * 2 ) > bzahl ) {
		c.hi++;
		c.low = bzahl - c.low;
		r_signum *= -1;
		}
	c.hi *= q_signum;
	c.low *= r_signum;
	return(c);
}

/*}}}  */
/*{{{  euclid, static*/
/*
*
* static elt_div_pair euclid(a,b);
*
*                                 /a\   /ggt(a,b)\
* berechnet Matrix A so, dass A * | | = |        |
*                                 \b/   \    0   /
*
 */
static elt_div_pair 
euclid (  /* frei nach gap */
    int a,
    int b
)
{
elt_div_pair r;
elt_div_prod temp;
int s, t, ha, hb;

  temp.hi  = 0;
  temp.low = 0;
  t        = 0;
  r.f1     = 0;
  r.f2     = 0;
  r.g1     = 0;
  r.g2     = 0;
  s  = (a >= 0) ? 1 : -1;
  hb= abs(b);
  ha= abs(a);
  while (hb != 0) {
    temp = ZZ_mindivmod(ha,hb);
    int q        = temp.hi;
    temp.hi  = 0;
    ha       = hb;
    hb       = temp.low;
    temp.low = 0;
    int t1 = t;
    t  = -t*q + s;
    s  = t1;
  }
  r.g1 = t;
  r.f1 = s;
  s    = ha - s * a;
  r.f2 = s/b;
  t   *= a;
  r.g2 = - t/b;
  return r;
}

/*}}}  */

/*{{{  elt_func, static*/
static matrix_TYP *
elt_func (matrix_TYP *mat)
{  
int min_sum, **E, h,*v, _min_;
int min_col, min_row, rM, cM, i, j, k,l;
boolean rc_min, flag;
#ifdef DEBUG
char *temp = "tempA";
#endif
matrix_TYP *elt;
elt_div_pair paar;

  if ( !(mat->flags.Integral & 1)) {
    fprintf (stderr, "Elementary divisors of a rational matrix ??\n");
    exit (3);
  }
  /*------------------------------------------------------------*\
  | Copy mat->array.SZ into E                    |
  \*------------------------------------------------------------*/
  rM = mat->rows;
  cM = mat->cols;
  min_sum = 0;
  elt = copy_mat(mat);
  E = elt->array.SZ;
  rc_min = TRUE;
  for ( i = 0; i < cM; i++) {
#ifdef DEBUG
    if (FALSE)
      put_mat(elt,temp,NULL,0);
    temp[4]++;
#endif
    /*---------------------------------------------------------*\
    |  Extended Version: Try to decrease the entries of the     |
    |  matrix by using linear combinations of rows and columns  |
    \*---------------------------------------------------------*/
    if (E[i][i] != 0) {
      for(j = i+1; j < rM; j++) {
        if((h= min_div(E[j][i],E[i][i])) != 0) {
          if( h == 1 ) {
            for(k = i; k < cM; k++) {
              if (E[i][k] != 0) {
                E[j][k] -= E[i][k];
              }
            }
          } else if (h == -1) {
            for(k = i; k < cM; k++) {
              if (E[i][k] != 0) {
                E[j][k] += E[i][k];
              }
            }
          } else {
            for(k = i; k < cM; k++) {
              if (E[i][k] != 0) {
                E[j][k] -= h*E[i][k];
              }
            } 
          }
        }
      }
      for(j = i+1; j < cM; j++) {
        if((h= min_div(E[i][j],E[i][i])) != 0) {
          if(h == 1) {
            for(k = i; k < rM; k++) {
              if (E[k][i] != 0) {
                E[k][j] -= E[k][i];
              }
            }
          } else if (h == -1) {
            for(k = i; k < rM; k++) {
              if (E[k][i] != 0) {
                E[k][j] += E[k][i];
              }
            } 
          } else {
            for(k = i; k < rM; k++) {
              if (E[k][i] != 0) {
                E[k][j] -= h*E[k][i];
              }
            }
          }
        }
      }
    }
    /* put_mat(elt,NULL,NULL,0); */
  /*{{{  */
    /*----------------------------------------------------------*\
    | Second extension: Trying to decrease the Entries of the    |
    | matrix by using a modificated pair Reduction.              |
    | (But it seems to make it slower in some cases)
    \*----------------------------------------------------------
    if ((i%3) == 0) {
      for(l = i+1; (l < rM) && (l < cM); l++) {
        if(E[l][l] != 0) {
      for(j = i+1; j < rM; j++)
        if(j != l) {
        if((h= min_div(E[j][l],E[l][l])) != 0) {
          if(h == 1)
          for(k = i; k < cM; k++) {
            if (E[l][k] != 0) E[j][k] -= E[l][k];
            }
          else if (h == MINUS_1)
          for(k = i; k < cM; k++) {
            if (E[l][k] != 0) E[j][k] += E[l][k];
            }
          else
          for(k = i; k < cM; k++) {
            if (E[l][k] != 0) E[j][k] -= h*E[l][k];
            }
          }
        }
      for(j = i+1; j < cM; j++)
        if(j != l) {
        if((h=  min_div(E[l][j],E[l][l])) != 0) {
          if(h == 1)
          for(k = i; k < rM; k++) {
            if (E[k][l] != 0) E[k][j] -= E[k][l];
            }
          else if (h == MINUS_1)
          for(k = i; k < rM; k++) {
            if (E[k][l] != 0) E[k][j] += E[k][l];
            }
          else
          for(k = i; k < rM; k++) {
            if (E[k][l] != 0) E[k][j] -= h*E[k][l];
            }
          }
          }
          }
        }
     put_mat(elt,NULL,NULL,2);
      }*/
  /*}}}  */
    /*---------------------------------------------------------*\
    | Find the minimum of all abs. values of non-zero entries  |
    | in the remaining (rM - i)-dimensional lower right      |
    | submatrix                           |
    \*---------------------------------------------------------*/
  
    min_row = min_col = i;
    _min_ = 0;
    for ( j = i;  (j < rM); j++ ) {
      for (k = i;  (k < cM); k++) {
        if ( E[j][k] != 0 ) {
          if (_min_ == 0) {
            _min_= abs(E[j][k]);
            min_row = j;
            min_col = k;
            if (rc_min) {
              min_sum = 0;
              for(l = i; l < rM; l++) {
                if (E[l][k] < 0) {
                  min_sum -= E[l][k];
                } else {
                  min_sum += E[l][k];
                }
              }
            }
          } else {
            flag = abs(E[j][k]) == abs(_min_) ? 0 : (abs(E[j][k]) < abs(_min_) ? -1 : 1 );
            if ( flag == -1 ) {
              _min_= abs(E[j][k]);
              if ((rc_min) && (k !=  min_col)) {
                min_sum = 0;
                for(l = i; l < rM; l++) {
                  if (E[l][k] < 0) {
                    min_sum -= E[l][k];
                  } else {
                    min_sum += E[l][k];
                  }
                }
              }
              min_row = j;
              min_col = k;
            } else if((rc_min) && (flag == 0) && (k!= min_col)) {
              int sum = 0;
              for(l = i; l < rM; l++) {
                if (E[l][k] < 0) {
                  sum -= E[l][k];
                } else {
                  sum += E[l][k];
                }
              }
              if(sum < min_sum) {
                min_row = j;
                min_col = k;
                min_sum = sum;
                sum = 0;
              }
            }
          }
        }
      }
    }
    if( _min_ == 0) {
      goto ende;     
    }
    if( _min_ == 1) {
      flag = FALSE; 
    } else {
      flag = TRUE;
    }
    rc_min = TRUE; /* (min & 0); */
  
    /*---------------------------------------------------------*\
    | Swap the min-entry into upper left corner            |
    \*---------------------------------------------------------*/
    v          = E[i];
    E[i]       = E[min_row];
    E[min_row] = v;
    for ( j = i; j < rM; j++ ) {
      h             = E[j][i];
      E[j][i]       = E[j][min_col];
      E[j][min_col] = h;
    }
  
    /*------------------------------------------------------*\
    | try to decrease the min using linear combinations    |
    | of rows                        |
    \*------------------------------------------------------*/
  marke:                                             
    if ( flag ) {
      for(j = rM-1; flag && (j > i); j--) {
        if ( (h = E[j][i]%_min_) != 0 ) {
          paar = euclid (E[i][i], E[j][i]);
          /* find differences for momo
          sprintf(comment,"%d %d", i,j);
          put_mat(elt,NULL,comment,0); */
  
          for ( k = i; k < cM ; k++) {
            h += E[i][k]*paar.f1+E[j][k]*paar.f2;
            E[j][k]= E[j][k]*paar.g2+E[i][k]*paar.g1;
            E[i][k] = h;
          }
          if ( (_min_= abs(E[i][i])) == 1 ) {
            flag = FALSE;                   
          }
        }
      }
    }
    if ( flag ) {
    /*------------------------------------------------------*\
    | try to decrease the min using linear combinations    |
    | of columns                      |
    \*------------------------------------------------------*/
      for ( j = cM-1; j > i; j--) {
        if ( (h= E[i][j]%_min_) != 0 ) {
          paar = euclid (E[i][i], E[i][j]);
          for ( k = i; k < rM; k++) {
            h      += E[k][i]*paar.f1+E[k][j]*paar.f2;
            E[k][j] = E[k][j]*paar.g2 + E[k][i]*paar.g1;
            E[k][i] = h;
          }
          if ( (_min_ = abs(E[i][i])) == 1 ) {
            flag = FALSE;                    
          }
        }
      }
    }
    /*---------------------------------------------------------*\
    | Clear the i-th column                     |
    \*---------------------------------------------------------*/
    for(j = i+1; j < rM ; j++ ) {
      h=E[j][i]/E[i][i];
      if (h != 0) {
        for ( k = i; k < cM ; k++ ) {
          if(E[i][k] != 0) {
            E[j][k] -= h*E[i][k];
          }
        }
      }
      if(E[j][i] != 0) {
        goto marke;        
      }
    }
    /*---------------------------------------------------------*\
    | Clear the i-th row                     |
    \*---------------------------------------------------------*/
    E[i][i] = abs(E[i][i]);
    for (j = i+1; j < cM; j++) {
      E[i][j] = 0;
    }
  }
ende:
  for (i = 1; i < elt->cols; i++) {
    if ( (E[i][i] % E[i-1][i-1]) != 0 ) {
      h= E[i-1][i-1];
      E[i-1][i-1]= GGT( E[i][i], h);
      h /= E[i-1][i-1];
      E[i][i] *= h;
      if (i > 1) {
        i -= 2;  
      }
    }
  }
  for (i = 0; i < elt->cols; i++) {
    if ( E[i][i] < 0 ) {
      E[i][i] = -E[i][i];
    }
  }
  
  for (i = 1; i < elt->cols; i++) {
    elt->array.SZ[0][i] = elt->array.SZ[i][i];
  }
  for (i = 1; i < elt->rows; i++) {
    free(elt->array.SZ[i]);
    if( elt->array.N  != NULL ) {
      free(elt->array.N[i]);
    }
  }
  elt->rows = 1;
  Check_mat(elt);
  return (elt);
}

/*}}}  */
/*{{{  elt_div*/
/*
@-----------------------------------------------------------------
@ matrix_TYP *elt_div(Mat)
@ matrix_TYP *Mat;
@
@ computes elementary divisors of Mat.
@ Result is a matrix with a single row that contains the
@ divisors.
@-----------------------------------------------------------------
*/
matrix_TYP *
elt_div (matrix_TYP *Mat)
{
matrix_TYP *Elt, *tmp;

  if(null_mat(Mat)) {
    Elt = init_mat(1,1,"0");
    return(Elt);
  }
  
  if (Mat->rows >= Mat->cols) {
    Elt = ggauss(Mat);
    tmp = tr_pose(Elt);
    free_mat(Elt);
    tgauss(tmp);
  } else {
    tmp = tr_pose(Mat);
    Elt = ggauss(tmp);
    free_mat(tmp);
    tmp = tr_pose(Elt);
    free_mat(Elt);
    tgauss(tmp);
  }
  Elt = elt_func(tmp);
  free_mat(tmp);
  return(Elt);
}
/*}}}  */
