#include "typedef.h"
#include "getput.h"
#include "matrix.h"
#include "dertype.h"
#include "tietzetrans.h"

extern int* multiply( matrix_TYP*, int*);
extern int* init_vector(int);

int is_equal_vec(vec_a, vec_b, dim)
int* vec_a;
int* vec_b;
{
int i;
for (i=0;i<dim;i++)
    {
    if (vec_a[i]!=vec_b[i]) return (0);
    }
return(1);
}

void split_relator(presentation, relatorno, beginsub, endsub)
presentation_TYP* presentation;
int relatorno;
int beginsub;
int endsub;
{
int i;
int k;
if (presentation->norelators>=presentation->ext_factor*EXT_SIZE)
   {
   presentation->ext_factor++;
   presentation->relators = (relator_TYP*) realloc (presentation->relators,
     sizeof(relator_TYP)*EXT_SIZE*presentation->ext_factor);
   }
k=presentation->norelators;
/* neuen Relator einfuegen */
presentation->relators[k].lhsnproduct = endsub - beginsub;
presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int)*(endsub-beginsub));
/* inserted tilman 21/3/97 */
presentation->relators[k].rhsproduct = NULL;
presentation->relators[k].rhsnproduct = 0;

for (i=beginsub;i<endsub;i++)
    {
    presentation->relators[k].lhsproduct[i-beginsub] = presentation->relators[relatorno].lhsproduct[i+1];
    }
/* alten Relator kuerzen */
for (i=endsub+1;i<presentation->relators[relatorno].lhsnproduct;i++)
    {
    presentation->relators[relatorno].lhsproduct[i-endsub+beginsub] =
    presentation->relators[relatorno].lhsproduct[i];
    }
/* hier evt. ein realloc einfuegen, falls es sich als notwendig erweist */
presentation->relators[relatorno].lhsnproduct = presentation->relators[relatorno].lhsnproduct - endsub + beginsub;
presentation->norelators++;
}

void sub_relator_search(presentation,gen_vector,dim)
presentation_TYP* presentation;
int* gen_vector;
int dim;
{

char name[1000];
matrix_TYP *tmp;

matrix_TYP **vector_list;
int i;
int j;
int finish;
int k;
int flag;

flag = 0;
for (i=(-1);i<presentation->norelators-1;)
    {
    if (!flag) 
       {
        i++;
       }
    else
       {
        flag =0;
       }
    vector_list = (matrix_TYP **) malloc (sizeof(matrix_TYP *)*presentation->relators[i].lhsnproduct);
    vector_list[0] = presentation->generators->list[
                     presentation->relators[i].lhsproduct[0]]->element;
    for (j=1;(j<presentation->relators[i].lhsnproduct)&&(!flag);j++)
        {
        vector_list[j] = mat_mul(vector_list[j-1],presentation->
          generators->list[presentation->relators[i].lhsproduct[j]]->element);
        for (k=0;((k<j)&&(!flag));k++)
            {
            if (mat_comp(vector_list[k], vector_list[j], dim)==0)
               {
             /* Subrelator gefunden ! */
              /* printf(" Subrelator gefunden, relator %d, start %d, ende %d \n",
                      i,k,j);
               */
             /* subrelator wird als neuer Relator eingefuegt, der alte
               relator wird entsprechend gekuerzt, der Durchlauf wird
               abgebrochen (durch setzen des flag) und derselbe Relator
               jetzt kuerzer wird erneut ueberprueft. */
               split_relator(presentation,i,k,j);
               flag = 1;
               finish = j;
               }
            }
        }
    if (!flag)
        {
        finish =presentation->relators[i].lhsnproduct;
        } 
    for (j=1;j<finish;j++)
        {
        free_mat(vector_list[j]);
        }
    free(vector_list);
    }
/* free(gen_vector); */
} 
