/******neukand berechnet kuerzeste Vektoren in vi+2T bzgl. gegebener Form A 
        form = neugram(mat,form);
       matrix_TYP *neukandidat_vectors(mat,form)   -> mat=T, form = TAT^tr ****/

/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: neukandidat.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include"typedef.h"
#include"matrix.h"
#include"presentation.h"

static int n;
static int laenge( int *vect, matrix_TYP *form);
static int zweihochk(int k);
static int *mach_restvec(int i);


/* berechne laenge von vec bzgl. des skalarprodukts form d.h. vec*form*vec^tr 
*/

static int laenge(vect,form)
int *vect;
matrix_TYP *form;

{	
	int i,j,l;
	l = 0;
 		for(i = 0; i < n; ++i)
		   for(j = 0; j < n; ++j)
			l += vect[j] * form->array.SZ[j][i] * vect[i];
		return l;
}

/*berechne 2^k */

static int zweihochk(k)
int k;
{
	int i, help;
 	help = 1;
	for(i = 0; i < k; ++i)
		help = help*2;
	return help;	
}

/*bilde Restklassenvertreter von T nach 2T*/

static int *mach_restvec(i)
int i;

{
     int k, rest;
     int *restvec;

/*Speicherplatz allocieren und restvec = (0,0,...,0) setzen*/

if((restvec = (int*)malloc(n * sizeof(int))) == 0)
	exit (2);
for (k= 0; k < n; ++k)
{ 
	restvec[k] = 0;
}

     rest = i;
	for(k = n-1; k > -1; --k)
	{
		if(zweihochk(k) <= rest)
			restvec[k] = 1;
		else
			restvec[k] = 0;
		rest -= restvec[k]*zweihochk(k);
	}
        if(rest != 0)
          printf("Fehler in mach_restvec mit argument i = %d\n", i);
	return(restvec);
}

/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *neukandidat_vectors(mat,Form)
@ matrix_TYP *mat,*Form;
@
@ calculates the shortest vectors in mat+T where T is the lattice determined by 
@ the form Form. The result is given in the rows of the output matrix.
@
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *neukandidat_vectors(mat,Form)
matrix_TYP *mat,*Form;
{
   int i,j,d;
   int *restvec;
   matrix_TYP *form;
   matrix_TYP *SV, *cand, *otto, *restvek;
   int candsize, esize;

   n = mat->rows;
   esize = 30;
   cand = init_mat(1, mat->cols+1, "");
   cand->array.SZ = (int **)realloc(cand->array.SZ, esize *sizeof(int *));
   candsize = esize;
   cand->rows = 0;

        form = neugram(mat,Form);
	sscal_mul(form,4);
	Check_mat(form);
      
   for(i = 1; i < zweihochk(mat->rows); ++i)
   {
	restvec = mach_restvec(i);
        d = laenge(restvec,form); 
	        restvek = init_mat(1,n,"");
		for(j=0;j<n;j++)
 		   restvek->array.SZ[0][j] = restvec[j];

		restvek->kgv = 2;  		/* da bzgl Basis von 2T geg. */
		sscal_mul(restvek,-1);  
	        Check_mat(restvek);

/*modshort liefert M={t in 2T | d(t-rv,0)<k} suche aber x aus rv+2T mit 
  d(x,0)=d(t+rv,0)<k, also modshort mit -restvek aufrufen, 
  ->M={t aus 2T | d(t+rv,0)<k}, x=rv+t, t aus M liefert d(x,0)=d(t+rv,0)<k  
*/

/*bestimme kuerzeste Vektoren in restvek+2T, starte mit max=d */

	SV = modshort_vectors(form, restvek,d);


/*habe zwei kuerzeste Vektoren in restvec+2T gefunden, falls(SV->rows == 1) und
  zwar SV->[0] und sein negatives 
  => dann ist SV relevant */

        if(SV->rows == 2)  
        {

	/* SV bzgl Basis von 2T geg, restvec ganzzahlig bzgl T, 
  	 d.h. 1/2*restvec ist Koordinatendarstellung bzgl. 2T ...
  	 bilde SV + vi */


	   for(j=0;j<SV->cols;j++)
           {
		SV->array.SZ[0][j]=2 * SV->array.SZ[0][j]; 
           	SV->array.SZ[0][j] += (SV->kgv) * restvec[j];
	   }
	   Check_mat(SV);


           if(cand->rows == 0)
           {
             for(j=0;j<=mat->rows;j++)
               cand->array.SZ[0][j] = SV->array.SZ[0][j];
              cand->rows++;
           }
           else
           {
              if(cand->rows == candsize)
              {
                 candsize += esize;
                 cand->array.SZ = 
                       (int **)realloc(cand->array.SZ, candsize *sizeof(int *));
              }
              cand->array.SZ[cand->rows] = 
                       (int *)malloc((mat->rows+1) *sizeof(int));

              for(j=0;j<=mat->rows;j++)
                cand->array.SZ[cand->rows][j] = SV->array.SZ[0][j];
              cand->rows++;
           }
        }/** if(SV->rows == 2)  **/
        free_mat(SV); SV = NULL;
        free_mat(restvek); restvek = NULL;
        free(restvec);
   }/** for(i) **/

   free_mat(form); form = NULL;
   cand->cols--;	/* letzte Spalte enthaelt Laenge */
   Check_mat(cand);

   otto = mat_mul(cand,mat); /*Ergebnis bzgl. Standardbasis*/
   Check_mat(otto);

   free_mat(cand); cand = NULL;

   return(otto);

}/**neukandidat_vectors()**/
