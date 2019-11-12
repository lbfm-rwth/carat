#include "typedef.h"
#include "utils.h"
#include "matrix.h"
#include "orbit.h"
#include "symm.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  vor_vertices.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP **voronoi_vertices(form, grp, anz, form_min, SV_no)
@ matrix_TYP *form;
@ bravais_TYP *grp;
@ int *anz, *form_min, *SV_no;
@
@ calculates the voronoi_vertices, i.e. the matrices in the space of
@ symmetric matrices, that are invariant under the action of G^{tr},
@ where the voronoi vertices are defined by the sum over all g in the
@ group 'grp'  of
@      g x^{tr}x g,
@ where x is a shortest vector of 'form'.
@
@ The number of these vertices is returned via (int *anz)
@ the number of shortest vectors of 'form' via (int *SV_no) and
@ the minimum of 'form' via (int *form_min)
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP **
voronoi_vertices (matrix_TYP *form, bravais_TYP *grp, int *anz, int *form_min, int *SV_no)
{
    int min_norm, orbit_anz;
    int *subdiv, *svi;
    matrix_TYP *SV, **vf;
    int i,j,k,l, n;

    n = form->cols;
    SV = shortest(form, &min_norm);
    *SV_no = SV->rows;
    *form_min = min_norm;
    subdiv = orbit_subdivision(SV, grp, &orbit_anz);
    vf = (matrix_TYP **)xmalloc(orbit_anz *sizeof(matrix_TYP *));
    for(i=0;i<orbit_anz;i++)
    {
      vf[i] = init_mat(n, n, "");
      vf[i]->flags.Symmetric = TRUE;
    }
    for(i=0;i<SV->rows;i++)
    {
      svi = SV->array.SZ[i];
      j = subdiv[i]-1;
      for(k=0;k<n;k++)
        for(l=0;l<=k;l++)
        {
           vf[j]->array.SZ[k][l] += (svi[k] * svi[l]);
           if(k!= l)
             vf[j]->array.SZ[l][k] = vf[j]->array.SZ[k][l];
        }
    }
    for(i=0;i<orbit_anz;i++)
    {
      l = 0;
      for(j=0;j<n && l != 1;j++)
         for(k=0;k<=j && l != 1 ;k++)
         {
           if(vf[i]->array.SZ[j][k] != 0)
           {
             if(l == 0)
               l = vf[i]->array.SZ[j][k];
             else
               l = GGT(l, vf[i]->array.SZ[j][k]);
             if(l<0)
              l = -l;
           }
         }
       if(l != 1)
       {
         for(j=0;j<n;j++)
          for(k=0;k<=j;k++)
          {
            vf[i]->array.SZ[j][k] /= l;
            vf[i]->array.SZ[k][j] = vf[i]->array.SZ[j][k];
          }
       }
    }
    free(subdiv);
    free_mat(SV);
    *anz = orbit_anz;
    return(vf);
}
