#include"typedef.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: init_bravais.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *init_bravais(dim)
@ int dim;
@
@ returns a bravais_TYP 'B' with B->dim = dim
@ and all other pointers an integers in B are set to NULL resp. 0
@ and B->divisors[i] = 0 for 0<= i < 100
@---------------------------------------------------------------------------
@
\**************************************************************************/


bravais_TYP *
init_bravais (int dim)
{
   int i;
   bravais_TYP *B;

   /* changed on 1/08/97 from malloc to calloc for paranoia reasons */
   if( (B = (bravais_TYP *)calloc(1,sizeof(bravais_TYP))) == NULL)
   {
     printf("malloc of 'B' in 'init_bravais' failed\n");
     exit(2);
   }
   B->dim = dim;
   for(i=0;i<100;i++)
     B->divisors[i] = 0;
   B->order = 0;
   B->gen_no = B->zentr_no = B->form_no = B->normal_no = B->cen_no = 0;
  B->gen = NULL;
  B->form = NULL;
  B->zentr = NULL;
  B->normal = NULL;
  B->cen = NULL;
  return(B);
}
