#include<typedef.h>
#include<getput.h>
#include<matrix.h>
#include<longtools.h>
#include<tools.h>
#include"zass.h"





matrix_TYP *scalar(long n,long a)
/* liefert eine skalarmatrix mit dem eintrag a */
{
  long i;
  matrix_TYP *erg;

  erg = init_mat(n,n,"k");

  for (i=0;i<n;i++){erg->array.SZ[i][i]=a;}

  return erg;
}   /* scalar */

matrix_TYP *matrizen_in_word(matrix_TYP **mat,matrix_TYP **matinv,word g)
/* geht davon aus, das alle (benoetigten) matrizen in **mat
   quadratisch von gleicher groesse sind. anderfalls wird nicht abgebrochen */
{
   matrix_TYP *erg,*tmp_matrix;
   long i,n;

   n = mat[0]->cols;
   erg = scalar(n,g.faktor);

   for (i=0;i<=g.last;i++){
     if (g.pointer[i]>0){
       tmp_matrix = mat_mul(erg,mat[g.pointer[i]-1]);
       free_mat(erg);
       erg = tmp_matrix;
     }
     else if (g.pointer[i]<0){
       /* 1ster fall, die inverse ist noch nicht berechnet */
       if (matinv[-g.pointer[i]-1]==NULL){
          matinv[-g.pointer[i]-1]= mat_inv(mat[-g.pointer[i]-1]);
       }
       tmp_matrix = mat_mul(erg,matinv[-g.pointer[i]-1]);
       free_mat(erg);
       erg = tmp_matrix;
     }
   }

   if (INFO_LEVEL & 4){
      put_mat(erg,NULL,NULL,2);
   }

   return erg;
}  /* matrizen_in_word */



static void norm_word(word *g)
/* soll die auftretenden worte vereinfachen */
{
  int i,j;

  /* wenn der faktor des relators = 0, dann auch der relator */
  if (g[0].faktor==0){
    g[0].last = -1;
    for (i=0;i<g[0].speicher;i++){
      g[0].pointer[i]=0;
    }
  }

  /* suchen nach einsen (representiert durch nullen) im relator */
  i = 0;
  while(i<=g[0].last){
    if (g[0].pointer[i]==0){
      for (j=i+1;j<=g[0].last;j++){
	g[0].pointer[j-1]=g[0].pointer[j];
      }
      g[0].last--;
    }
    else{
       i++;
    }
  }

  /* suchen nach ausdruecken wie x_j * x_j^(-1) */
  i=0;
  while(i<g[0].last){
    if (g[0].pointer[i]==(-g[0].pointer[i+1])){
      g[0].pointer[i]=0;
      g[0].pointer[i+1]=0;
      /* und jetzt die enstandenen nullen eliminieren
	 (geht nicht durch vertauschen der while-schleifen, weil
	 so x_1 * x_2 * x_2^(-1) * x_1^(-1) nicht erkannt wuerde) */
      norm_word(g);
      /* nur der sicherheit halber */
      i= -1;
    }
    i++;
  }

  return;
}

int wordfree(word *a)
{

  free(a[0].pointer);

  return TRUE;
}

static int ini_word(word *a,int laenge)
{
  a[0].pointer = (long *) malloc(laenge * sizeof(long));
  a[0].speicher = laenge;
  a[0].last = 0;

  return TRUE;
}


static int test_relator(matrix_TYP **mat,matrix_TYP **matinv,word *relator,long anz)
/* testet, ob die relationen in *relator erfuellt sind:
   wenn ja: return true; sonst false */
{
   long i,n;
   int flag;
   matrix_TYP *id,*tmp;

   flag = TRUE;
   i = 0;
   n = mat[0]->cols;
   id = scalar(n,1);

   while(i<anz && flag==TRUE){
     tmp = matrizen_in_word(mat,matinv,relator[i]);

     if (INFO_LEVEL & 4){
       put_mat(tmp,NULL,"tmp in test_relator",2);
       put_mat(id,NULL,"id in test_relator",2);
     }

     if (cmp_mat(tmp,id)!=0){
       flag = FALSE;
     }
     i++;
     free_mat(tmp);
   }

   free_mat(id);
   return flag;
}   /* test_relator */

void matrix_2_word(matrix_TYP *matrix,word *relator,long zeile)
{
  long i,hilf;

  if (relator[0].speicher !=0 ){
     free(relator[0].pointer);
     relator[0].speicher = 0;
  }

  ini_word(relator,matrix->cols);

  relator[0].last = matrix->cols -1;
  relator[0].faktor = 1;

  hilf = 0;
  for (i=0;i<=relator[0].last;i++){
    relator[0].pointer[i] = matrix->array.SZ[zeile][i+hilf];
    if (relator[0].pointer[i] == 0){i--;hilf++;relator[0].last--;}
  }

  norm_word(relator);

  return;
}  /* matrix_2_word */

static matrix_TYP *fox_deriv_mat(matrix_TYP **mat, matrix_TYP **matinv,word a,long j)
/* soll rekursiv die j-te fox-derivation vom word a berechnen und sie
   in erg abspeichern */
{
  word g,h;
  long i,teil;
  matrix_TYP *erg,
             *nu_g,
             *nu_h,
             *g_mat,
             *tmp;

  if (INFO_LEVEL & 4){
     printf("fox_deriv_mat\n");
  }

  norm_word(&a);

  if ( a.last==0 || a.last == (-1) ){
    if (a.pointer[0]==j){
       erg = scalar(mat[0]->cols,1);
    }
    else if (a.pointer[0]==(-j)){
       tmp = scalar(mat[0]->cols,-1);
       if (matinv[j-1]==NULL){
          matinv[j-1]=mat_inv(mat[j-1]);
       }
       erg = mat_mul(tmp,matinv[j-1]);
       free_mat(tmp);
    }
    else{
       erg = scalar(mat[0]->cols,0);
    }
  }
  else{

  /* die ersten teil erzeuger von a werden in g
     gespeichert, der rest in h  */
    teil = 1;
    ini_word(&g,teil);
    ini_word(&h,a.last);

    for (i=0;i<=a.last;i++){
      if (i<teil){
        g.pointer[i] = a.pointer[i];
      }
      else{
	h.pointer[i-teil] = a.pointer[i];
      }
    }

    g.last= teil-1;
    h.last= a.last-teil;
    g.faktor=1;
    h.faktor=1;
    nu_g = fox_deriv_mat(mat,matinv,g,j);
    nu_h = fox_deriv_mat(mat,matinv,h,j);

    g_mat = matrizen_in_word(mat,matinv,g);

    tmp = mat_mul(g_mat,nu_h);

    erg = mat_add(nu_g,tmp,One,One);

    wordfree(&g);
    wordfree(&h);

    free_mat(g_mat);
    free_mat(nu_g);
    free_mat(nu_h);
    free_mat(tmp);

  }

  if (INFO_LEVEL &4){
     printf("return aus fox_deriv_mat\n");
  }

  return erg;
}


matrix_TYP *calc_B(matrix_TYP **mat,long anz_erzeuger)
{
  long i,k,l,n;
  matrix_TYP *B;

  n = mat[0]->cols;

  /* reservieren des speichers fuer die matrix B, die die
     erzeuger - id uebereinander enthaelt */
  B = init_mat(n*anz_erzeuger,n,"k");

  /* belegen der matrix B */
  for (i=0;i<anz_erzeuger;i++){
    for (k=0;k<n;k++){
      for (l=0;l<n;l++){
        if (k == l){
          B->array.SZ[k+i*n][l]=mat[i]->array.SZ[k][l]-1;
        }
        else{
          B->array.SZ[k+i*n][l]=mat[i]->array.SZ[k][l];
        }
      }
    }
  }

  if (INFO_LEVEL & 4){
    put_mat(B,NULL,"%B",2);
  }

  return B;
}   /* calc_B */

static matrix_TYP *calc_A(matrix_TYP **mat,matrix_TYP **matinv,word *relator,
                                       int erzeuger,int relatoren)
{
  matrix_TYP *eingesetzt,
             *A;
  long i,j,k,l,n;

  /* reservieren von speicher fuer die grosse matrix A */
  n = mat[0]->cols;
  A = init_mat(n*relatoren,n*erzeuger,"k");

  /* berechnen der j-ten foxableitung des i-ten relators */
  for (i=0;i<relatoren;i++){
    for (j=1;j<=erzeuger;j++){

      eingesetzt = fox_deriv_mat(mat,matinv,relator[i],j);

      if (INFO_LEVEL & 4){
         put_mat(eingesetzt,NULL,NULL,2);
      }

      Check_mat(eingesetzt);

      /* umspeichern der eintraege aus der fox-ableitung in das
         gleichungssystem A */
      for (k=0;k<n;k++){
        for (l=0;l<n;l++){
          A->array.SZ[n*i+k][n*(j-1)+l] = eingesetzt->array.SZ[k][l];
        }
      }

      /* die matrix eingesetzt wird wieder benutzt */
      free_mat(eingesetzt);
    }
  }

  if (is_option('s')){
    put_mat(A,FILENAMES[2],"system of congruences to determine H^1(G,**)",2);
  }

  if (INFO_LEVEL & 4){
    put_mat(A,NULL,NULL,2);
  }

  return A;
}   /* calc_A  */



matrix_TYP** cohomology(long *dim,
                        matrix_TYP **mat,
                        matrix_TYP **matinv,
                        word *relator,
                        int erzeuger,
                        int relatoren)
{   matrix_TYP *A,
               *B,
               *B_tr,
               *cozykel,
               *elementar,
               **tmp;

    int corang_a,rang_b,erg,i;

    if (INFO_LEVEL &4){
       printf("in cohomology\n");
    }

    /* ueberprueft die gegeben relatoren */
    if (test_relator(mat,matinv,relator,relatoren)==FALSE){
       fprintf(stderr,"Error in Cohomlogy:\n");
       fprintf(stderr,"One of the given relators it not satisfied\n");
       fprintf(stderr,"by this generating set.\n");
       fflush(stderr);
       exit(3);
    }

    /* berechnet die matrix des fuer die cozykel zu
       loesenden gleichungssystems  */
    A = calc_A(mat,matinv,relator,erzeuger,relatoren);

    if (INFO_LEVEL & 4){
       put_mat(A,NULL,"A",2);
    }

    /* enthaelt nur die erzeuger - id hintereinander,
       spannt also den raum der coraender auf */
    B = calc_B(mat, erzeuger);
    B_tr = tr_pose(B);


    /* shrink the linear system of equations by a gauss algorithm,
       note that this also removes the dependence of the result
       from the given presentation */
    long_row_hnf(A);

    /* loesen der gleichungssysteme */
    tmp = cong_solve(A);

    cozykel = tmp[0];
    elementar = tmp[1];

    if (INFO_LEVEL &4){
       put_mat(cozykel,NULL,"cozykel",2);
       put_mat(elementar,NULL,"elementar",2);
    }

    rang_b = tgauss(B_tr);

    /* berechnen des coranges von A */
    corang_a = 0;
    for (i=0;i<elementar->cols;i++){
      if (elementar->array.SZ[i][i]==0){
         corang_a++;
      }
    }

    dim[0] = erg = corang_a - rang_b;

    if (INFO_LEVEL & 4){
       printf("corang_a %d\n",corang_a);
       printf("rang_b %d\n",rang_b);
       printf("erg %d\n",erg);
    }
 
    if (erg == 0){
       /* die Q-dimensionen der Coraender und Cozykel sind gleich, d.h.
          es kommt nur torsion dazu */
       for (i=elementar->cols-1;0<=i;i--){
         if (elementar->array.SZ[i][i]==0){
            kill_col(elementar,i);
            kill_row(elementar,i);
            kill_col(cozykel,i);
         }
       }
    }
    else{
       /* es interesiert */
       for (i=elementar->cols-1;0<=i;i--){
         if (elementar->array.SZ[i][i]!=0){
            kill_col(elementar,i);
            kill_row(elementar,i);
            kill_col(cozykel,i);
         }
       }
    }

    if (INFO_LEVEL & 4){
       put_mat(A,NULL,NULL,2);
       put_mat(B,NULL,NULL,2);
       put_mat(elementar,NULL,NULL,2);
       put_mat(cozykel,NULL,NULL,2);
    }

    free_mat(A);
    free_mat(B);
    free_mat(B_tr);
    free_mat(elementar);

    /* free_mat(cozykel); */
    /* free_mat(tmp2); */
    tmp[0] = cozykel;
    tmp[1] = tmp[3];

    tmp = (matrix_TYP **) realloc(tmp,3 *sizeof(matrix_TYP *));
    /* tmp[1] = elementar; */

    return tmp;
}  /* cohomology  */

