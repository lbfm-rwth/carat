#include "typedef.h"
#include "matrix.h"
#include "getput.h"
#include "datei.h"
#include "orbit.h"
#include "longtools.h"
#include "sort.h"
#include "bravais.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: gittstabneu.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *gittstab(grp, X)
@ bravais_TYP *grp;
@ matrix_TYP *X;
@
@   gittstab calculates the stabilizer of a lattice X under the        |
@   operation of the group 'grp' acting via left-multiplication        |
@---------------------------------------------------------------------------
@
\**************************************************************************/
bravais_TYP *gittstabneu(grp, X)
bravais_TYP *grp;
matrix_TYP *X;
{
   int i,
       j,
       search,
       length,
       options[6];

   bravais_TYP *stab;

   matrix_TYP **orbit,
               *INV,
               *Y;

   /* make the options for orbit */
   for (i=0;i<6;i++) options[i] = 0;
   options[1] = 4;
   options[3] = 1;

   /* gauss reduce a copy of X */
   Y = copy_mat(X);
   long_col_hnf(Y);

   /* reserve memory as asked for in orbit_alg */
   stab = init_bravais(grp->dim);

   orbit = orbit_alg(Y,grp,stab,options,&length);

   /* throw away the orbit */
   for (i=0;i<length;i++){
      free_mat(orbit[i]);
   }
   free(orbit);
   free_mat(Y);

   /* try to reduce the number of generators */
   /* firstly sort them */
   mat_quicksort(stab->gen,0,stab->gen_no-1,mat_comp);
   for (i=0;i<stab->gen_no-1;i++){
      if (mat_comp(stab->gen[i],stab->gen[i+1]) == 0){
         free_mat(stab->gen[i+1]);
         stab->gen_no--;
         for (j=i+1;j<stab->gen_no;j++){
            stab->gen[j] = stab->gen[j+1];
         }
         i--;
      }
   }

   /* second pass, use the inverses */
   /* bare in mind that the generator list is sorted and does not have
      duplicates */
   for (i=0;i<stab->gen_no;i++){
      INV = mat_inv(stab->gen[i]);
      search = mat_search(INV,stab->gen,stab->gen_no,mat_comp);
      free_mat(INV);
      /* watch for elements of order 2 !!! */
      if (search != -1 && search != i){
         if (search<i){
            fprintf(stderr,"Error in gittstabneu\n");
            exit(3);
         }
         free_mat(stab->gen[search]);
         stab->gen_no--;
         for (j=search;j<stab->gen_no;j++){
            stab->gen[j] = stab->gen[j+1];
         }
      }
   }

   /* reallocate the memory for stag->gen */
   stab->gen = (matrix_TYP **) realloc(stab->gen,
                                       sizeof(matrix_TYP *)*stab->gen_no);

   /* calculate the order of the stabilizer (if possible) */
   if (grp->order != 0){
      stab->order = grp->order/length;
      factorize_new(stab->order,stab->divisors);
   }
   else{
      stab->order = 0;
      for (i=0;i<100;i++) stab->divisors[i] = 0;
   }

   return(stab);
}
