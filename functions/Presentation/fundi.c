#include "typedef.h"
#include "polyeder.h"
#include "dertype.h"
#include "roundtype.h"

/* Prototypenliste ---------------------------------------------------------*/
int get_smallest (wall_list_TYP*);
fund_domain* fundamental_domain( wall_list_TYP*, derivedsg_TYP* );
wall_TYP** convert_walls( wall_list_TYP* );
/*---------------------------------------------------------------------------*/
extern fund_domain* first_fuber(wall_TYP**, int);
extern void free_wall (wall_TYP **);
extern void put_fund_domain(fund_domain*);
extern int refine_fund_domain(fund_domain*, wall_TYP*);
extern void free_mat(matrix_TYP*);
extern matrix_TYP* copy_mat( matrix_TYP* );
extern matrix_TYP* mat_mul( matrix_TYP*, matrix_TYP* ); 
extern matrix_TYP* init_mat (int,int,char*);
extern matrix_TYP* mat_inv( matrix_TYP*);

/*-----------------------------------------------------------------------------
   Funktion, die aus der Liste der Waende die Wand mit dem kleinsten Abstand
   zum generating vector heraussucht. Da es sich um einen sortierten Baum
   handelt (sortiert nach eben diesen Abstaenden) braucht dabei nicht die
   gesamte Liste abgelaufen werden. 
   Parameter:
       walls - Liste der Waende.
   Rueckgabewert:
       Index der gefundenen Wand in der Liste.
-----------------------------------------------------------------------------*/
int get_smallest (walls)
wall_list_TYP* walls;
{
int lauf=0, laufalt=0, laufuralt, mark=(-1), markalt=0;
while (lauf!=nichts)
   {
   laufuralt = laufalt;
   laufalt = lauf;
   if (walls->wall_list[lauf]->right == nichts)
      {
      markalt = laufuralt;
      mark = lauf;
      lauf = nichts;
      }
   else
      lauf = walls->wall_list[lauf]->right;
   }
if (mark==(-1))
   {
   walls->wall_list[laufuralt]->right=nichts;
   }
else
   {
   walls->wall_list[markalt]->right = walls->wall_list[mark]->left; 
   walls->wall_list[mark]->left = nichts;
   laufalt=mark;
   }
return (laufalt);
}

/*-----------------------------------------------------------------------------
    Konvertierung der Walls ist noetig, da refine_fuber und first_fuber ein
    etwas anderes Format benutzen.
    Parameter:
         walls - Liste der Waende im Format wall_Typ (definiert in roundtype.h)
    Rueckgabewert:
         Liste der Waende im Format wall_TYP (definiert in 
             ...kristall/include/typedef.h)
-----------------------------------------------------------------------------*/
wall_TYP** convert_walls( walls )
wall_list_TYP* walls;
{
wall_TYP** allwalls;
int i,j,isml;
allwalls = (wall_TYP**) malloc (sizeof(wall_TYP*) * walls->firstfree);
/* Triviale Wand (Nullvektor) rausschmeissen */
isml = get_smallest (walls);
for (i=1;i<walls->firstfree;i++)
    {
    allwalls[i-1] = (wall_TYP*) malloc (sizeof(wall_TYP));
    allwalls[i-1]->gl = (int*) malloc (sizeof(int)*walls->dimension);
    isml = get_smallest (walls);
    allwalls[i-1]->nproduct = walls->wall_list[isml]->nproduct;
    allwalls[i-1]->product = (int*) malloc (sizeof(int)*allwalls[i-1]->nproduct);
    for (j=0;j<allwalls[i-1]->nproduct;j++)
        allwalls[i-1]->product[j] = walls->wall_list[isml]->product[j];
    for (j=0;j<walls->dimension;j++)
        {
        allwalls[i-1]->gl[j] = walls->wall_list[isml]->hplane[j];
        }
    allwalls[i-1]->norm = walls->wall_list[isml]->dist;
    allwalls[i-1]->dim = walls->dimension;
    allwalls[i-1]->mat = NULL;
    }
return(allwalls);
}

/*-----------------------------------------------------------------------------
   Bestimmen eines Fundamentalbereichs.
   Parameter:
              walls - Waende um einen gegebenen generating vector
              generators - Generatoren der Gruppe
   Rueckgabewert: 
              Der Fundamentalbereich.
-----------------------------------------------------------------------------*/
fund_domain* fundamental_domain( walls, generators )
wall_list_TYP* walls;
derivedsg_TYP* generators;
{
wall_TYP** allwalls;
int wallno;
int i, j;
int erg;
fund_domain* fudo;
matrix_TYP* hilf;
allwalls = convert_walls( walls);
wallno = walls->firstfree-1;
/* fudo = first_fuber (allwalls, wallno); */
fudo = first_polyeder (allwalls, wallno);

if (fudo == 0) 
	{
	printf(" Error: first_fuber couldn' t build a fundamental domain !\n\n");
	exit(3);
	}
fprintf(stderr," First fundamental domain has %d walls. \n",fudo->wall_no);
for (i=0;i<wallno;i++)
    {
    /* don't bore the user */
    if (i % 100 == 0){
       fprintf(stderr," . ");
    }

    /* if (refine_fund_domain(fudo, allwalls[i])==1) */
    if (refine_polyeder(fudo, allwalls[i])==1)
       allwalls[i]=NULL;
    }

fprintf(stderr,"\n");
fprintf(stderr," Refined fundamental domain has %d walls. \n",fudo->wall_no);
/* zu den Waenden gehoerende Matrizen erzeugen */
for (i=0;i<fudo->wall_no;i++)
    {
    fudo->wall[i]->mat = init_mat(walls->dimension, walls->dimension,"1");

    for (j=0;j<fudo->wall[i]->nproduct;j++)
        {
        hilf = mat_mul(generators->list[fudo->wall[i]->product[j]]->element,
                       fudo->wall[i]->mat);
        free_mat(fudo->wall[i]->mat);
        fudo->wall[i]->mat=hilf;
        }
    /*fudo->wall[i]->mat = mat_inv(hilf);*/

    /* deleted the next two lines because they are useless and cause trouble
       if fudo->wall[i]->nproduct == 0 ,  tilman 6/3/97
    fudo->wall[i]->mat = copy_mat(hilf);
    free_mat(hilf); */
    }
/* nicht mehr benoetigten Speicherplatz freigeben  */
for (i=0;i<wallno;i++)
    {
    if (allwalls[i]!=NULL)
        {
        if (allwalls[i]->gl!=NULL) free(allwalls[i]->gl);
        if (allwalls[i]->nproduct!=0) free(allwalls[i]->product);
        free(allwalls[i]);
        }
    }
/*put_fund_domain(fudo);*/
free(allwalls);
return(fudo);
}
