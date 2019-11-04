#include "typedef.h"
#include "matrix.h"
#include"getput.h"
#include"datei.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: free_bravais.
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_bravais(grp)
@ bravais_TYP *grp;
@
@ frees the space allocated for the pointers of the bravais_TYP *grp.
@---------------------------------------------------------------------------
@
\**************************************************************************/

void 
free_bravais (bravais_TYP *grp)
{
int i;

  if ( grp ) {
    if ( grp->gen_no != 0 &&  grp->gen ) {
      for(i=0; i< grp->gen_no; i++) {
          free_mat(grp->gen[i]);
      }  
      free ( grp->gen );
    }
    
    if ( grp->form_no != 0 && grp->form ) {
      for(i=0; i< grp->form_no; i++) {
          free_mat(grp->form[i]);
      }
      free( (int *)grp->form);
    }
    
    if ( grp->zentr_no != 0  && grp->zentr ) {
      for(i=0; i< grp->zentr_no; i++)
      {
          free_mat(grp->zentr[i]);
      }
      free( (int *)grp->zentr);
    }
      
    if ( grp->normal_no != 0  && grp->normal ) {
      for(i=0; i< grp->normal_no; i++)
      {
          free_mat(grp->normal[i]);
      }
      free( (int *)grp->normal);
    }
    
    if ( grp->cen_no != 0  && grp->cen ) {
      for(i=0; i< grp->cen_no; i++)
      {
        free_mat(grp->cen[i]);
      }
      free( (int *)grp->cen );
    }

    free( (int *)grp);
  }

}


/*{{{}}}*/
