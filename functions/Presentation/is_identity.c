#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"

extern int* copy_product_r (int*, int);
extern matrix_TYP* mat_inv (matrix_TYP*);
extern matrix_TYP* copy_mat (matrix_TYP*);
extern matrix_TYP* mat_mul (matrix_TYP*, matrix_TYP*);
extern derivedsg_TYP* init_derivedgen (int);
extern void free_mat(matrix_TYP*);
extern derived_TYP* copy_derived(derived_TYP*);

int is_identity (presentation, dim)
presentation_TYP* presentation;
int dim;
{
static matrix_TYP* mat;
matrix_TYP* del;
int i,j;
if (presentation->relators[presentation->norelators].lhsnproduct==1)
   {
   mat = copy_mat(presentation->generators->list[
             presentation->relators[presentation->norelators].lhsproduct[0]
                              ]->element);
   return(0);
   }
else
   {
   del = mat;
   mat = mat_mul (del, presentation->generators->list[
        presentation->relators[presentation->norelators].lhsproduct[
        presentation->relators[presentation->norelators].lhsnproduct-1]
                ]->element);
   free_mat(del);
   for (i=0;i<dim;i++)
       {
       for (j=0;j<dim;j++)
           {
           if (mat->array.SZ[i][j]!=delta(i,j)) 
              {
              return(0);
              }
           }
       }
   free_mat(mat);
   return(1);
   }
}

int is_identity_new (presentation,anzahl, dim)
/* Beachte: mat behaelt den Wert vom vorherigen Aufruf. Sollte
   diese Funktion beim letzten Aufruf nicht 1 zurueckgeben, so bleibt 
   der Speicherplatz von mat belegt. */
presentation_TYP* presentation;
int anzahl;
int dim;
{
static matrix_TYP* mat;
matrix_TYP* del;
int i,j;
if (presentation->relators[presentation->norelators].lhsnproduct==anzahl)
   {
   mat = copy_mat(presentation->generators->list[
             presentation->relators[presentation->norelators].lhsproduct[0]
                              ]->element);
   for (i=0;i<(anzahl-1);i++)
       {
       del = mat;
/* es kann ruhig die Matrix des letzten Elements im Produkt multipliziert
   werden, denn die letzten anzahl-Elemente sind gleich. */
       mat = mat_mul (del, presentation->generators->list[
        presentation->relators[presentation->norelators].lhsproduct[
        /* presentation->relators[presentation->norelators].lhsnproduct-1] */
        i+1]
                ]->element);
       free_mat(del);
       } 
   return(0);
   }
else
   {
   for (i=0;i<anzahl;i++)
       {
       del = mat;
       mat = mat_mul (del, presentation->generators->list[
        presentation->relators[presentation->norelators].lhsproduct[
        presentation->relators[presentation->norelators].lhsnproduct-1]
                ]->element);
       free_mat(del);
       }
   for (i=0;i<dim;i++)
      {
      for (j=0;j<dim;j++)
          {
          if (mat->array.SZ[i][j]!=delta(i,j)) 
             {
             return(0);
             }
          }
      }
  free_mat(mat);
  return(1);
  }
}
