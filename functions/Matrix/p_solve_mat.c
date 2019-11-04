#include "typedef.h"
#include "tools.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: p_solve.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP **p_solve (anz, L_mat, R_mat, option)
@ int *anz;
@ matrix_TYP **L_mat, **R_mat ;
@ int option;
@
@ Solves L_mat[i] * X = X * R_mat[i] modulo act_prime simultaneous
@  for 0 <= i<= anz.
@ The return is a basis of the space of solution.
@ the dimension of this space is returned by the pointer anz.
@ 
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  p_solve*/
matrix_TYP **
p_solve (int *anz, matrix_TYP **L_mat, matrix_TYP **R_mat, int option)
{
/*{{{  variables*/
int **M, **N, *v, f, help, flag ;
int *M_j, *M_act, numax;
matrix_TYP **E;
boolean L_null, R_null ;
int *ss, *nuin;
int i, j, jj, k, kk, l, cR, cL, cM, rR, rL, rM, rN, rg = 0, act, p;
int **TS;

/*}}}  */
  /*{{{  var init & error-check*/

  p = act_prime;
  L_null = (L_mat == NULL);
  R_null = (R_mat == NULL);

  cL = L_null ? 1 : L_mat[0]->cols;
  rL = L_null ? 1 : L_mat[0]->rows;
  rR = R_null ? 1 : R_mat[0]->rows;
  cR = R_null ? 1 : R_mat[0]->cols;

  if(!(L_null || R_null)) {
    if((cL != rL) || (rR != cR)) {
      fprintf(stderr,"Error in input: Matrizes of wrong dimension\n");
      exit(3);
    }
  }

  cM = cL * rR;
  rM = *anz*rL*cR;
  rN = rL * cR;
  rg = *anz*cR;    

  /*}}}  */
  /*{{{  alloc memory*/
  N    = (int **)malloc(rN*sizeof(int *));
  M    = (int **)malloc(rM*sizeof(int *));
  nuin = (int *) malloc((cM+1)*sizeof(int));
  ss   = (int *) malloc(cM*sizeof(int));

  /*}}}  */
  
  for ( i = 0  ; i < *anz; i++) {
  /*{{{  */
    if ( !L_null ) {
      if ( (L_mat[i]->cols != cL) || (L_mat[i]->rows != rL)) {
        fprintf (stderr, "Error in input\n");
        exit (3);
      }
    }
    if ( !R_null ) {
      if (( R_mat[i]->rows != rR ) || ( R_mat[i]->cols != cR )) {
        fprintf (stderr, "Error in input\n");
        exit (3);
      }
    }

    for ( j = 0; j < rN; j++ ) {
      N[j] = (int *)calloc(cM,sizeof(int));
    }

    if ( !L_null ) {
      /*{{{  */
      flag = L_mat[i]->prime != 0;
      TS = L_mat[i]->array.SZ;
      for ( j = 0, jj = 0; j < rL; j++, jj += rR) {
        for ( k = 0, kk = 0; k < cL; k++, jj -= rR) {
          if(flag) {
            help = TS[j][k];
          } else if ((help = TS[j][k] % p) < 0) {
            help += p;
          }
          if(help) {
            for( l = 0; l < rR; l++,jj++, kk++) {
              N[jj][kk] = help;                 
            }
          } else {
            jj += rR;
            kk += rR;
          }
        } 
      }
      /*}}}  */
    }
    if ( !R_null ) {
      /*{{{  */
      flag = R_mat[i]->prime != 0;
      TS = R_mat[i]->array.SZ;
      for ( j = 0; j < cR; j++) {
        for ( k = 0; k < rR; k++) {
          if(flag) {
            help = TS[k][j];
          } else if ((help = TS[k][j] % p) < 0) {
            help += p;
          }
          if(help) {
            for(l=0,jj=j, kk=k; l<cL; l++, jj+=cR, kk+=rR) {
              N[jj][kk] = S(N[jj][kk],-help);
            }
          }
        }
      }
      /*}}}  */
    }

    /* permutierte Eingabe in Liste */
    for(j = 0; j < rL; j++) {
      for(k = 0; k < cR; k++) {
        M[j*rg+i*cR+k] = N[(rL-1-j)*cR+k];
      }
    }
  /*}}}  */
  }
  free(N);
  rg = 0;
    /*------------------------------------------------------------*\
    | Gauss Algorithm            |
    \*------------------------------------------------------------*/
  for ( i = 0; i < cM; i++ ) {
#if 0
  /*{{{  auskommentiert*/
    for(pri = 0; pri < rM; pri++) {
      for(prj = 0; prj < cM; prj++)
        if(M[pri][prj] != 0)
          printf("%d ",M[pri][prj]);
        else
          printf("  ");
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
  
  if ( (*anz = cM - rg) != 0 ) {
  /*{{{  */
    /*========================================================*\
    || Matrix has not full rank: clear it upwards             ||
    \*=========================================================*/
    for (i = rg-1; i > 0; i--) {
      /*{{{  */
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
      /*}}}  */
    }
    if((option & 2) || (!L_null && !R_null) ) {
      /*{{{  */
      E = (matrix_TYP **)malloc(*anz*sizeof(matrix_TYP *));
      for ( i = 0 , k = 0, l = 0; i < cM; i++ ) {
        if ( i == ss[k] ) {
          k ++;
          continue;
        }
        E[l] = init_mat(cL, rR, "ip");
        E[l]->prime = act_prime;
        N = E[l]->array.SZ;
        for ( j = 0; j < k; j ++) {
          if ( M[j][i] != 0 ) {
            N[ss[j]/rR][ss[j]%rR] = p-M[j][i];
          }
        }
        N[i/rR][i%rR] = 1;
        l ++;
      }
      /*}}}  */
    } else {
      /*{{{  */
      E = (matrix_TYP **)malloc(1*sizeof(matrix_TYP *));
      E[0] = init_mat(*anz, cM,"ip");
      E[0]->prime = act_prime;
      N = E[0]->array.SZ;
      for ( i = 0 , k = 0, l = 0; i < cM; i++ ) {
        if ( i == ss[k] ) {
          k ++;
          continue;
        }
        for ( j = 0; j < k; j ++) {
          if ( M[j][i] != 0 ) {
            N[l][ss[j]] = p-M[j][i];
          }
        }
        N[l][i] = 1;
        l ++;
      }
      *anz = 1;
      /*}}}  */
    }
  /*}}}  */
  } else {
    E = NULL;
  }
  
  /*{{{  free memory*/
  for (i = 0; i < rM; i++ ) {
    free(M[i]);
  }
  free(M);
  free(ss);
  free(nuin);
  /*}}}  */
  
  return (E);
}
/*}}}  */
