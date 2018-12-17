#include "typedef.h"
#include "longtools.h"
#include "matrix.h"
#include "sort.h"
#include "getput.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: orbit_alg.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP **orbit_alg(M, G, S, option, length)
@ bravais_TYP *G, *S;
@ matrix_TYP *M;
@ int *option, *length;
@ 
@ 'orbit_alg calculates the orbit of the matrix M under the group G.
@ 
@  The length of the orbit is returned by the pointer 'length'.
@  It is possible to calculate the stabiliser of M in G.
@  It is returned in the pointer 'S' which has to be allocated as '*bravais_TYP'
@  before calling the function.
@ 
@  'option' is an array of size 5.
@  The following options are possible:
@ 
@  option[0] = 0: the group operates from the left: x->gx
@  option[0] = 1: the group operates from the right: x->xg
@  option[0] = 2: the group operates via x-> g x g^(tr) (from the left)
@  option[0] = 3: the group operates via x-> g^(tr) x g (from the right)
@  option[0] = 4: the group operates by conjugation: x-> g x g^(-1)
@                 (from the left)
@  option[0] = 5: the group operates by conjugation: x-> g^(-1) x g
@                 (from the right)
@ 
@  option[1] = 1: the group operates on the pairs {M, -M} 
@  option[1] = 2: the group  caluluates the orbit of the set of rows of M;
@  option[1] = 3: combination of 1 and 2.
@  option[1] = 4: the group operates on sublattices of Z^n.
@
@  option[2] = 0: the whole orbit is calculated
@  option[2] = n: only the first n elements of the orbit are calculated
@
@  option[3] = 1: the stabiliser is calculated.
@
@  option[4] = 0: the full stabiliser is calculated.
@  option[4] = n: only the first n elements of the stabiliser are calculated
@  option[4] = -n: only the first n elements of the stabiliser are calculated,
@                  afterwards the algorithm is stopped.
@
@  option[5] = 1: Also the inverse of the generators are used in the algorithmn.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*****************************************
   orbit_alg bestimmt die Bahn der Matrix M unter der Gruppe G

   Die Bahnlaenge wird ueber den Zeiger length zurueckgegeben.
   Soll der Stabilisator von M in G berechnet werden, so muessen
   die Inversen der Erzeuger von G in Ginv angegeben werden.
   In diesem Fall wird der Stabilisator ueber den Zeiger S zurueckgegeben.
   Der Zeiger S muss aber bereits als *bravais_TYP allokiert sein.

   option ist ein array der Groesse 5.
   Folgende Optionen sind moeglich:

   option[0] = 0: die Gruppe operiert durch Linksmultiplikation: x->gx
   option[0] = 1: die Gruppe operiert durch Rechtsmultiplikation: x->xg
   option[0] = 2: die Gruppe operiert durch x-> g x g^(tr) (von links)
   option[0] = 3: die Gruppe operiert durch x-> g^{tr} x g (von rechts)
   option[0] = 4: die Gruppe operiert durch Konjugationation: x-> g x g^(-1)
                  (von links)
   option[0] = 5: die Gruppe operiert durch Konjugationation: x-> g^(-1) x g
                  (von rechts)

   option[1] = 1: die Gruppe operiert auf den Paaren {M, -M} 
   option[1] = 2: die Gruppe  berechnet die Bahn der Menge der Zeilen von M;
   option[1] = 3: Kombination von 1 und 2.
   option[1] = 4: die Gruppe operiert auf Teilgittern von Z^n.

   option[2] = 0: die volle Bahn wird berechnet
   option[2] = n: die ersten n Elemente der Bahn werden berechnet

   option[3] = 1: der Stabilisator wird berechnet

   option[4] = 0: der volle Stabilisator wird berechnet
   option[4] = n: die ersten n Elemente des Stabilisators werden berechnet
   option[4] = -n: die ersten n Elemente des Stabilisators werden berechnet,
                   dannach wird der Algorthmus abgebrochen.

   option[5] = 1: Auch die Inversen der Erzeuger werden im Algorithmus
                  angewendet.
******************************************/


static matrix_TYP **Erz;
static matrix_TYP **Erz_inv;
static int *orbit_opt;
static matrix_TYP *hash_mat;
static int hash_prime = 130003;


static void matrix_speichern(mat, L, listsize, anz)
matrix_TYP *mat;
matrix_TYP ***L;
int *listsize, *anz;
{
 int i, j, k;

extern matrix_TYP *init_mat();
   if((*anz) == 0)
      if ((*L=(matrix_TYP **)malloc((*listsize)*sizeof(matrix_TYP *)))==NULL)
      {
         fprintf (stderr, "malloc failed\n");
         exit (2);
      }
   if ((*anz) >= (*listsize))
   {
     *listsize += EXT_SIZE;
     if ((*L=(matrix_TYP **)realloc(*L,(*listsize)*sizeof(matrix_TYP *)))==NULL)
     {
        fprintf (stderr, "realloc failed\n");
        exit (2);
     }
   }
   (*L)[(*anz)] = mat;
   (*anz) ++;
}



static struct baum *addbaum(mat, L, anz, verz, schonda)
matrix_TYP *mat;
matrix_TYP **L;
int anz;
struct baum *verz;
int *schonda;
{
  int vergleich;
  *schonda = -1;
  if(verz == NULL)
  {
    verz = (struct baum *) malloc(sizeof(struct baum));
    verz->no = anz;
    verz->left = verz->right = NULL;
  }
  else
  {
     vergleich = mat_comp(mat, L[verz->no]);
     if(vergleich<0)
       verz->left = addbaum(mat, L, anz, verz->left, schonda);
     if(vergleich > 0)
       verz->right = addbaum(mat, L, anz, verz->right, schonda);
     if(vergleich == 0)
      *schonda = verz->no;
  }
  return(verz);
}


struct baum *hash_addbaum(mat, L, anz, verz, schonda, hashnumber, hashverz)
matrix_TYP *mat;
matrix_TYP **L;
int anz;
struct baum *verz;
int *schonda, hashnumber, *hashverz;
{
  int vergleich;
  *schonda = -1;
  if(verz == NULL)
  {
    verz = (struct baum *) malloc(sizeof(struct baum));
    verz->no = anz;
    verz->left = verz->right = NULL;
  }
  else
  {
     if(hashnumber != hashverz[verz->no])
     {
        if(hashnumber < hashverz[verz->no])
         verz->left = hash_addbaum(mat, L, anz, verz->left, schonda, hashnumber, hashverz);
        else
         verz->right = hash_addbaum(mat, L, anz, verz->right, schonda, hashnumber, hashverz);
     }
     else
     {
       vergleich = mat_comp(mat, L[verz->no]);
       if(vergleich<0)
         verz->left = hash_addbaum(mat, L, anz, verz->left, schonda, hashnumber, hashverz);
       if(vergleich > 0)
         verz->right = hash_addbaum(mat, L, anz, verz->right, schonda, hashnumber, hashverz);
       if(vergleich == 0)
        *schonda = verz->no;
     }
  }
  return(verz);
}




int *make_orbit_options()
{
   int *option;

   option = (int *)malloc(6 *sizeof(int));
   option[0] = 0;
   if(is_option('l') == TRUE)
       option[0] = 0;
   if(is_option('r') == TRUE)
       option[0] = 1;
   if(is_option('t') == TRUE)
       option[0] = 2;
   if(is_option('t') == TRUE && is_option('r') == TRUE)
       option[0] = 3;
   if(is_option('k') == TRUE)
       option[0] = 4;
   if(is_option('k') == TRUE && is_option('r') == TRUE)
       option[0] = 5;

   option[1] = 0;
   if(is_option('p') == TRUE)
      option[1] = 1;
   if(is_option('u') == TRUE)
      option[1] = 2;
   if(is_option('p') == TRUE && is_option('u') == TRUE)
      option[1] = 3;
   if(is_option('g') == TRUE)
      option[1] = 4;

   option[2] = optionnumber('L');
   option[3] = is_option('S');
   option[4] = optionnumber('S');
   if (option[4] == -1) option[4] = 0;
   option[5] = is_option('i');
 return(option);
}

static void qswap (mat, k, l)
matrix_TYP *mat;
int k,l;
{  int *temp;

   temp            = mat->array.SZ[k];
   mat->array.SZ[k] = mat->array.SZ[l];
   mat->array.SZ[l] = temp;
   temp = NULL;
}




static void orbit_qsort (mat, left, right, orbit_comp)
matrix_TYP *mat;
int left, right;
int (*orbit_comp)();
{  
   int i, last;
   void qswap ();

   if ( left >= right )
      return;
   qswap (mat, left, (left + right) / 2);
   last = left;
   for ( i  = left+1;
         i <= right ;
         i ++        )
      if ( (*orbit_comp)(mat, i, left) < 0 )
         qswap (mat, ++last, i);
   qswap (mat, left, last);
   orbit_qsort (mat, left, last - 1, orbit_comp);
   orbit_qsort (mat, last + 1, right, orbit_comp);
}

static int orbit_comp (mat, k, l)
matrix_TYP *mat;
int k,l;
{  int i;
   int **M;
   int cM;

   M  = mat->array.SZ;
   cM = mat->cols;
      for ( i  = 0 ; i<cM; i++)
      {
         if( M[k][i] > M[l][i])
               return (-1);
         else
           if( M[k][i] < M[l][i])
                  return (1);
      }
      return (0);
}



static void standartisieren(mat)
matrix_TYP *mat;
{
  int   i,
        j,
        k;

  if(orbit_opt[1] == 3)
  {
    for(i=0; i<mat->rows; i++)
    {
       k=0;
       while(k<mat->cols && mat->array.SZ[i][k] == 0)
          k++;
       if(k<mat->cols && mat->array.SZ[i][k] <0)
       {
          for(j=k; j<mat->cols; j++)
             mat->array.SZ[i][j]  = -mat->array.SZ[i][j];
       }
    }
  }
  if(orbit_opt[1] == 1)
  {
    k = 0; i = 0;
    while(i<mat->rows && k == 0)
    {
      j=0;
      while(j<mat->cols && k==0)
      {
          if(mat->array.SZ[i][j] < 0)
            k= -1;
          if(mat->array.SZ[i][j] > 0)
            k= 1;
          j++;
      }
      i++;
    } 
    if(k == -1)
    {
      for(i=0;i<mat->rows;i++)
        for(j=0;j<mat->cols;j++)
         mat->array.SZ[i][j]  = -mat->array.SZ[i][j];
    }
  }
  if(orbit_opt[1] == 2 || orbit_opt[1] == 3)
    orbit_qsort(mat, 0, mat->rows-1, orbit_comp);
  if(orbit_opt[1] == 4)
  {
     if(orbit_opt[0] == 0 || orbit_opt[0] == 2 || orbit_opt[0] == 4){
       long_col_hnf(mat);
     }
     if(orbit_opt[0] == 1 || orbit_opt[0] == 3 || orbit_opt[0] == 5)
       long_row_hnf(mat);
  }
}


static matrix_TYP *tr_mul(A, B)
matrix_TYP *A, *B;
{
  int i, j, k;
  matrix_TYP *erg;

  extern matrix_TYP *init_mat();
  erg = init_mat(A->rows, B->rows, "");
  for(i=0; i<erg->rows; i++)
  {
    for(j=0; j<erg->cols; j++)
    {
      erg->array.SZ[i][j] = 0;
      for(k=0; k<A->cols; k++)
        erg->array.SZ[i][j] += A->array.SZ[i][k] * B->array.SZ[j][k];
    }
  }
  erg->kgv = A->kgv *B->kgv;
  return(erg);
}

static matrix_TYP *grp_mul(A,B)
matrix_TYP *A, *B;
{
  matrix_TYP *erg;
  extern matrix_TYP *mat_mul();
  if(orbit_opt[0] == 0 || orbit_opt[0] == 2 || orbit_opt[0] == 4)
    erg = mat_mul(B,A);
  if(orbit_opt[0] == 1 || orbit_opt[0] == 3 || orbit_opt[0] == 5)
    erg = mat_mul(A,B);
  return(erg);
}

static matrix_TYP *operation(M, i)
matrix_TYP *M;
int i;
{
   matrix_TYP *erg, *waste, *waste1;
 
   extern matrix_TYP *mat_mul();
   extern matrix_TYP *tr_mul();

   if(orbit_opt[0] == 0)
      erg = mat_mul(Erz[i], M);
   if(orbit_opt[0] == 1)
      erg = mat_mul(M, Erz[i]);
   if(orbit_opt[0] == 2)
   {
      waste = tr_mul(M, Erz[i]);
      erg = mat_mul(Erz[i], waste);
      free_mat(waste);
   }
   if(orbit_opt[0] == 3)
   {
     waste = tr_pose(Erz[i]);
     waste1 = mat_mul(M, Erz[i]);
     erg = mat_mul(waste, waste1);
     free_mat(waste); free_mat(waste1);
   }
   if(orbit_opt[0] == 4)
   {
     waste = mat_mul(Erz[i], M);
     erg = mat_mul(waste ,Erz_inv[i]);
     free_mat(waste);
   }
   if(orbit_opt[0] == 5)
   {
     waste = mat_mul(Erz_inv[i], M);
     erg = mat_mul(waste ,Erz[i]);
     free_mat(waste);
   }
   return(erg);
}

static void make_hash_mat(M)
matrix_TYP *M;
{
   int i,j;
   extern matrix_TYP *init_mat();
  
   hash_mat = init_mat(M->rows, M->cols, "");
   for(i=0; i<M->rows;i++)
     for(j=0;j<M->cols;j++)
     {
         hash_mat->array.SZ[i][j] = rand();
         hash_mat->array.SZ[i][j] %= 1301;
     }
}


static int hash_number(M)
matrix_TYP *M;
{
    int i,j,h;
    h = 0;
    for(i=0;i<M->rows;i++)
      for(j=0;j<M->cols;j++)
        h = (h + M->array.SZ[i][j] * hash_mat->array.SZ[i][j])%hash_prime;
    return(h);
}

extern void free_baum(struct baum *p)
{

   if (p!= NULL){
      if (p->left != NULL){
         free_baum(p->left);
      }
      if (p->right != NULL){
         free_baum(p->right);
      }
      free(p);
   }


   return;
}

matrix_TYP **orbit_alg(M, G, S, option, length)
bravais_TYP *G, *S;
matrix_TYP *M;
int *option, *length;
{
  int i,j,k;
  matrix_TYP **erg;
  matrix_TYP **hist, **histinv;
  matrix_TYP *I, *A, *K1, *K2;
  int *hashnumbers;
  int h, *hashverz;
  int ergsize, histsize, histinvsize, Ssize;
  int erganz, histanz, histinvanz;
  int schonda, sch, abbruch;
  int num, erzj, erz_anz;
  int sopt;
  struct baum *Sverz;
  struct baum *ergverz;

  extern matrix_TYP *init_mat();
  extern matrix_TYP *mat_inv();
  extern matrix_TYP *einheitsmatrix();
  extern void matrix_speichern();

  ergsize = 0;
  histsize = 0;
  histinvsize = 0;
  erganz = 0;
  histanz = 0;
  histinvanz = 0;
  Ssize = 0;
  ergverz = NULL;
  Sverz = NULL;
  orbit_opt = option;
  erzj = 0;
  if(orbit_opt[5] == 1)
    erz_anz = 2 * G->gen_no;
  else
    erz_anz = G->gen_no;
  Erz = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
  for(i=0;i< G->gen_no;i++)
    Erz[i] = G->gen[i];
  if(orbit_opt[5] == 1)
  {
    for(i=0;i<G->gen_no;i++)
      Erz[G->gen_no + i] = mat_inv(G->gen[i]);
  }
  if(orbit_opt[3] == TRUE || orbit_opt[0] == 4 || orbit_opt[0] == 5)
  {
    Erz_inv = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
    for(i=0;i<erz_anz;i++)
       Erz_inv[i] = mat_inv(Erz[i]);
  }
  make_hash_mat(M);


  standartisieren(M);
  h = hash_number(M);
  ergverz = hash_addbaum(M, erg, erganz, ergverz, &schonda, h, NULL);
  erg = (matrix_TYP **)malloc(EXT_SIZE *sizeof(matrix_TYP));
  erg[0] = init_mat(M->rows, M->cols, "");
  hashverz = (int *)malloc(EXT_SIZE *sizeof(int));
    hashverz[0] = h;
  for(i=0; i<M->rows; i++)
    for(j=0; j<M->cols; j++)
      erg[0]->array.SZ[i][j] = M->array.SZ[i][j];
  erg[0]->kgv = M->kgv;
  erganz = 1;
  ergsize = 100;
  if(orbit_opt[3] == 1)
  {
     I = einheitsmatrix(Erz[0]->cols);
     matrix_speichern(I, &hist, &histsize, &histanz);
     matrix_speichern(I, &histinv, &histinvsize, &histinvanz);
     S->gen_no = 0;
     S->dim = Erz[0]->cols;
  }
  sopt = orbit_opt[3];

  abbruch = FALSE;
  num = 0;
  while(abbruch == FALSE)
  {
     for(erzj=0; erzj<erz_anz && abbruch == FALSE; erzj++)
     {
        A = operation(erg[num], erzj);
        standartisieren(A);
        h = hash_number(A);
        ergverz = hash_addbaum(A, erg, erganz, ergverz, &schonda, h, hashverz);
        if(schonda != -1 && sopt == TRUE)
        {
           K1 = grp_mul(hist[num], Erz[erzj]);
           K2 = grp_mul(K1, histinv[schonda]);
           if(mat_comp(K2, I) != 0)
           {
            Sverz = addbaum(K2, S->gen, S->gen_no, Sverz, &sch);
            if(sch == -1)
             matrix_speichern(K2, &S->gen, &Ssize, &S->gen_no);
           }
           else{
              /* inserted this case on 21 Oct 98 tilman */
              free_mat(K2);
              K2 = NULL;
           }
           if ((schonda == -1 || (schonda != -1 && sch != -1)) && K2 != NULL)
              free_mat(K2);
           free_mat(K1); K2 = NULL;
        }
        if(schonda == -1)
        {
          if(erganz == ergsize)
          {
            if((hashverz = (int *)realloc(hashverz, (ergsize+EXT_SIZE)*sizeof(int))) == NULL)
            {
               printf("realloc failed\n");
               exit(2);
            }
          }
          hashverz[erganz] = h;
          matrix_speichern(A, &erg, &ergsize, &erganz);
          A = NULL;
          if(sopt == TRUE)
          {
             K1 = grp_mul(hist[num], Erz[erzj]);
             matrix_speichern(K1, &hist, &histsize, &histanz);
             K1 = NULL;
             K1 = grp_mul(Erz_inv[erzj], histinv[num]);
             matrix_speichern(K1, &histinv, &histinvsize, &histinvanz);
             K1 = NULL;
          }
        }
        if(schonda != -1)
          free_mat(A);
        if(orbit_opt[2] > 0 && erganz == orbit_opt[2])
          abbruch = TRUE;
        if(orbit_opt[4] <0 && -orbit_opt[4] == S->gen_no)
          abbruch = TRUE;
        if(orbit_opt[4] > 0 && orbit_opt[4] == S->gen_no)
         sopt = FALSE;
     }
     num++;
     if(num == erganz)
       abbruch = TRUE;
  }
  *length = erganz;
  if(orbit_opt[3] == TRUE)
  {
    for(i=1; i<histanz; i++)
    {
      free_mat(hist[i]);
      free_mat(histinv[i]);
    }
    free(hist); free(histinv);
    free_mat(I);
  }

  /* changed tilman 11/3/97 from
  if(orbit_opt[3] == TRUE || orbit_opt[0] == 3)
  to : */
  if(orbit_opt[3] == TRUE || orbit_opt[0] == 4)
  {
    for(i=0;i<erz_anz;i++)
       free_mat(Erz_inv[i]);
    free(Erz_inv);
  }
  if(orbit_opt[5] == 1)
  {
    for(i=0;i<G->gen_no;i++)
      free_mat(Erz[G->gen_no + i]);
  }
  free(Erz);
  free_mat(hash_mat);

  /* inserted 16/1/07 tilman */
  free(hashverz);
  free_baum(ergverz);
  free_baum(Sverz);

  /* inserted 16/07/97 tilman */
  for (i=0;i<length[0];i++)
    Check_mat(erg[i]);

  return(erg);
}

