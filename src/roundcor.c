/*
   Programm zur Berechnung der Mittelsenkrechten zwischen einem Vektor
   (durch die Funktion set_vector gesetzt; muss in allgemeiner Lage
    sein, d.h. Stabilisator ist das Einselement) und den Bahnelementen
   (der Bahn, die von diesem Vektor unter Gruppenoperation gebildet
    wird).
    Erstellt am: 05.10.94
    Autor: Joerg Kock, Lehrstuhl B Mathematik, RWTH-Aachen
    Aenderung am: 11.10.1994
    Zweck: Sortierung der Waende nach Abstand zum Vektor, um den der
           Dirichletbereich angelegt wurde.
    Aenderung am: 21.10.1994
    Zweck: Einbau des Aufrufs zur Bestimmung des Fundamentalbereichs
    (i.e. der Funktion "fundamental_domain", die mit Hilfe der Funktionen 
    first_fuber und refine_fund_domain einen Fundamentalbereich bestimmt
    (Funktionen aus: /usb/kristall/functions/fuber )).
    Aenderung am: 02.01.1995
    Zweck: Einbau des Aufrufs zur Bestimmung der Relationen "generate_pres"
           Einbau des Aufrufs "calc_orders" zur Bestimmung der Ordnungen
           der Generatoren.
    Aenderung am: 05.01.1995
    Zweck: Einbau des Aufrufs "tietzetrans". D.h. Bestimmung der Praesentation
    in den alten Generatoren.
    Aenderung zwischen 05.01.1995 und 01.03.1995
    Zweck: Verkuerzung der Relatoren mittels Subrelatorsearch. Bringt 
    Kuerzung der Gesamtlaenge bis zu 50% (im Durchschnitt etwa 15%).
*/

#include "typedef.h"
#include "bravais.h"
#include "getput.h"


#include "dertype.h"
#include "tietzetrans.h"
#include "roundtype.h"

extern int SCOUNT;
extern int INFO_LEVEL;
extern int SFLAG;

extern void old_as_new( derivedsg_TYP*, derivedsg_TYP*, int);
extern derivedsg_TYP* copy_derivedsg( derivedsg_TYP*);
extern void calc_orders (presentation_TYP*,int*, int);
extern void dumprelators();
extern matrix_TYP* fget_mat(FILE*);
extern matrix_TYP* init_mat(int, int, char*);
extern matrix_TYP* copy_mat(matrix_TYP*);
extern matrix_TYP* mat_mul(matrix_TYP*, matrix_TYP*);
extern bravais_TYP* get_bravais(char*);
extern derivedsg_TYP* init_derivedgen();
extern void dumpmat();
extern fund_domain* fundamental_domain(wall_list_TYP*, derivedsg_TYP*);
extern presentation_TYP* generate_pres(fund_domain*,int*);
extern matrix_TYP* mat_inv(matrix_TYP*);

/* Prototypenliste: */
void 		calc_walls (derivedsg_TYP*, int*, wall_list_TYP*, matrix_TYP*);
wall_Typ* 	wall_operation (derived_TYP*, int, wall_Typ*, matrix_TYP*, int* );
void 		add_vec( int*, matrix_TYP*);
int* 		multiply(matrix_TYP*, int*);
int* 		subtract_vec(int*, matrix_TYP*);
int* 		copy_product_r(int*, int);
wall_list_TYP* 	init_wall_list(int);
int 		is_same_wall(wall_Typ*, wall_Typ*, int);
int 		is_greater_wall(wall_Typ*, wall_Typ*);
int 		is_in_wall_list(wall_list_TYP*, wall_Typ* , int*);
void 		insert_wall(wall_list_TYP*, wall_Typ*);
void 		free_wall_l(wall_Typ*);
void 		free_wall_list(wall_list_TYP*);
void 		free_vector(int*);
int*		init_vector(int);
void 		set_vector(int*, int,int);
void 		convert_bravais(bravais_TYP*, derivedsg_TYP*);
matrix_TYP* 	mat_vec_mul( matrix_TYP*, int*);
void 		dump_walls(wall_list_TYP*);
matrix_TYP*	get_grammatrix(char*, int);
int		distance (int*, int*, int);
void 		free_derivedsgcomp (derivedsg_TYP*);
void 		dumpproduct (int*, int);


void 
dumpproduct (int *prod, int nprod)
{
int i;
for (i=0;i<nprod;i++)
    {
    printf(" %d ",prod[i]);
    }
printf("\n");
}

void 
free_presentation (presentation_TYP *presentation)
{
int i;
free_derivedsgcomp (presentation->generators);
for (i=0;i<presentation->norelators;i++)
    {
    if (presentation->relators[i].rhsnproduct!=0)
       {
       free(presentation->relators[i].rhsproduct);
       presentation->relators[i].rhsnproduct = 0;
       presentation->relators[i].rhsproduct = NULL;
       }
    if (presentation->relators[i].lhsnproduct!=0)
       {
       free(presentation->relators[i].lhsproduct);
       presentation->relators[i].lhsnproduct = 0;
       presentation->relators[i].lhsproduct = NULL;
       }
    }
free(presentation->relators);
free(presentation);
}

/*-------------------------------------------------------------------------
   Funktion, die eine Waendeliste erzeugt. 
   Parameter: Dimension des Raumes
--------------------------------------------------------------------------*/
wall_list_TYP *
init_wall_list (int dim)
{
wall_list_TYP* hilf;
hilf= (wall_list_TYP*) malloc (sizeof(wall_list_TYP));
hilf->wall_list=(wall_Typ**) malloc (EXT_SIZE*sizeof(wall_Typ*));
hilf->sizemult=1;
hilf->firstfree=0;
hilf->dimension=dim;
return(hilf);
}

/*-------------------------------------------------------------------------
   Funktion, die zwei Waende vergleicht und bei Gleichheit eine Eins zurueck-
   gibt, sonst Null.
   Parameter: Die zwei Waende die verglichen werden sollen und die Dimension
   des Raumes.   
--------------------------------------------------------------------------*/
int 
is_same_wall (wall_Typ *wall1, wall_Typ *wall2, int dim)
{
int i;
for (i=0;i<dim;i++)
    if (wall1->hplane[i]!=wall2->hplane[i]) return(0);
return(1);
}

/*-------------------------------------------------------------------------
   Alphabetisches" vergleichen zweier Waende. Gibt Eins zurueck, wenn
   wall1 groesser ist als wall2, sonst Null.
   Parameter: wall1 und wall2, die zu vergleichenden Waende. dim die 
   Dimension des Raumes. 
--------------------------------------------------------------------------*/
/* ausgetauscht gegen Funktion, die den Abstand vergleicht
int is_greater_wall(wall1, wall2, dim)
wall_Typ* wall1;
wall_Typ* wall2;
int dim;
{
int i;
for (i=0;i<dim;i++)
    if (wall1->hplane[i]>>wall2->hplane[i]) return(1);
return(0);
}
*/
int 
is_greater_wall (wall_Typ *wall1, wall_Typ *wall2)
{
int i;
if (wall1->dist>wall2->dist) return(1);
return(0);
}

/*-------------------------------------------------------------------------
   
--------------------------------------------------------------------------*/
/*int is_in_wall_list(walls, wall, index)
wall_list_TYP* walls;
wall_Typ* wall;
int* index;
{
int lauf=0;
while (lauf!=nichts)
     {
     if (is_same_wall(walls->wall_list[lauf], wall, walls->dimension))
        return(1);
     lauf = walls->wall_list[lauf]->next;
     }
return(0);
} auskommentiert, da lineare Liste zu langsam*/

/*-------------------------------------------------------------------------
    Funktion, die ueberprueft, ob die Wand wall bereits in der Liste walls
    vorhanden ist.   
--------------------------------------------------------------------------*/
int 
is_in_wall_list (wall_list_TYP *walls, wall_Typ *wall, int *index)
{
int lauf=0;
while (lauf!=nichts)
     {
     if (is_same_wall(walls->wall_list[lauf], wall, walls->dimension))
        return(1);
     if (walls->wall_list[lauf]->dist>wall->dist)
        lauf = walls->wall_list[lauf]->right;
     else
        lauf = walls->wall_list[lauf]->left;
     }
return(0);
}

/*--------------------------------------------------------------------------- 
   Funktion, die die Waende in einer Baumstruktur ablegt (zum schnelleren
   finden). Vergleichsoperation ist der Abstand zum Generierungsvektor, dessen
   wert in ->norm gespeichert ist. 
-----------------------------------------------------------------------------*/
void 
insert_wall (wall_list_TYP *walls, wall_Typ *newwall)
{
int lauf=0, laufalt;
/* Neuen Speicherplatz anfordern falls noetig */
if (walls->firstfree>=EXT_SIZE*walls->sizemult)
   {
   walls->sizemult++;
   walls->wall_list=(wall_Typ**) realloc(walls->wall_list,EXT_SIZE*sizeof(wall_Typ*)*walls->sizemult);
   }
if (walls->firstfree!=0)
   {
   do
      {
      laufalt = lauf;
      if (walls->wall_list[lauf]->dist>newwall->dist)
         lauf= walls->wall_list[lauf]->right;
      else
         lauf= walls->wall_list[lauf]->left;
      } while (lauf!=nichts);
   if (walls->wall_list[laufalt]->dist>newwall->dist)
      walls->wall_list[laufalt]->right = walls->firstfree;
   else
      walls->wall_list[laufalt]->left = walls->firstfree;
   }
newwall->left=nichts;
newwall->right=nichts;
walls->wall_list[walls->firstfree]=newwall;
walls->firstfree++;
}

/*-------------------------------------------------------------------------
   Einfuegen einer Wand in die Wandliste (lineare aufsteigend nach Abstand
   zum generierungsvektor sortierte Liste). 
--------------------------------------------------------------------------*/
/*void insert_wall(walls, newwall)
wall_list_TYP* walls;
wall_Typ* newwall;
{
int lauf=0, laufalt=0;
if (walls->firstfree>=EXT_SIZE*walls->sizemult)
   {
   walls->sizemult++;
   walls->wall_list=(wall_Typ**) realloc(walls->wall_list,EXT_SIZE*sizeof(wall_Typ*)*walls->sizemult);
   }
if (walls->firstfree!=0)
   {
   lauf = walls->first;
   while (lauf!=nichts)
       {
       if (is_greater_wall(walls->wall_list[lauf], newwall))
          {
          newwall->next = lauf;
          if (lauf!=walls->first)
              walls->wall_list[laufalt]->next = walls->firstfree;
          else
              walls->first = walls->firstfree;
          laufalt=nichts;
          lauf=nichts;
          }
       else
          {
          laufalt = lauf;
          lauf = walls->wall_list[lauf]->next;
          }
       }
    if (laufalt!=nichts)
       {
       newwall->next=nichts;
       walls->wall_list[laufalt]->next = walls->firstfree; 
       }
   }
else
   {
   walls->first=0;
   newwall->next=nichts;
   }
walls->wall_list[walls->firstfree]=newwall;
walls->firstfree++;
}
auskommentiert, weil lineare Liste zu langsam*/

/*-------------------------------------------------------------------------
    Freigeben des Speicherplatzes einer Wand 
--------------------------------------------------------------------------*/
void 
free_wall_l (wall_Typ *wall)
{
free(wall->hplane);
if (wall->nproduct!=0){
   free(wall->product);
   wall->product = NULL;
   wall->nproduct = 0;
}
free(wall);
}

/*-------------------------------------------------------------------------
   Freigeben des Speicherplatzes einer Waendeliste. 
--------------------------------------------------------------------------*/
void 
free_wall_list (wall_list_TYP *walls)
{
int i;
for (i=0;i<walls->firstfree;i++)
    {
    free_wall_l(walls->wall_list[i]);
    }
free(walls->wall_list);
free(walls);
}

/*-------------------------------------------------------------------------
    Vektor initialisieren. 
--------------------------------------------------------------------------*/
int *
init_vector (int dim)
{
int* hilf;
hilf= (int*)calloc (dim,sizeof(int));
return(hilf);
}

/*-------------------------------------------------------------------------
   Speicherplatz eines Vektors freigeben. 
--------------------------------------------------------------------------*/
void 
free_vector (int *vec)
{
free(vec);
}

/*-------------------------------------------------------------------------
   generating vector setzen. Es muss ein Vektor gefunden werden dessen
   Stabilisator Eins ist. Z.Z. wird ein beliebiger Vektor genommen und
   gehofft, dass er diese Eigenschaft hat. Ist der Stab nicht Eins bricht
   das Programm ab und es muss ein anderer Vektor ausprobiert werden. 
--------------------------------------------------------------------------*/
void 
set_vector (int *vec, int dim, int tries)
{
int i;

for (i=0;i<dim;i++)
   vec[i] = (2* (rand()%2) - 1) * (rand() % 5);

return;

vec[0]=2;
vec[1]=(-1);
for (i=2;i<dim;i++)
    {
    vec[i]=(-vec[i-2])/abs(vec[i-2])*(abs(vec[i-2])+1);
    }
}

/*-------------------------------------------------------------------------
   Konvertieren von bravais_TYP -> derivedsg_TYP 
--------------------------------------------------------------------------*/
void 
convert_bravais (bravais_TYP *brav, derivedsg_TYP *gen)
{
int i;
derived_TYP* hilf;
for (i=0;i<brav->gen_no;i++)
    {
    hilf=(derived_TYP*) malloc (sizeof(derived_TYP));
    hilf->element=copy_mat(brav->gen[i]);
    hilf->nproduct = 0;
    insert_list(gen,hilf);
    }
}

/*-------------------------------------------------------------------------
    Kopieren einer Integerliste 
--------------------------------------------------------------------------*/
int *
copy_product_r (int *prod, int nprod)
{
int* hilf;
int i;
hilf = (int*) malloc (sizeof(int)*(nprod+1));
for (i=0;i<nprod;i++)
    hilf[i] = prod[i];
return(hilf);
}

/*-------------------------------------------------------------------------
    Subtrahieren zweier Vektoren, die in unterschiedlichen Typen gespeichert
    sind. 
--------------------------------------------------------------------------*/
int *
subtract_vec (int *vec, matrix_TYP *matvec) 
{
int* hilf;
int i;
hilf = (int*) malloc (sizeof(int)*matvec->rows);
for (i=0;i<matvec->rows;i++)
   {
   hilf[i] = vec[i] - matvec->array.SZ[i][0];
   }
return(hilf);
}

/*-------------------------------------------------------------------------
    Multiplizieren einer Matrix mit einem Vektor (der nicht im matrix_TYP
    vorliegt). 
--------------------------------------------------------------------------*/
int *
multiply (matrix_TYP *mat, int *vec)
{
int* hilf;
int i, j, sum;
hilf = (int*) malloc (sizeof(int)*mat->rows);
for (i=0;i<mat->rows;i++)
    {
    sum = 0;
    for (j=0;j<mat->cols;j++)
        {
        sum+= mat->array.SZ[i][j]*vec[j]; 
        }
    hilf[i] = sum;
    }
return(hilf);
}

/*-------------------------------------------------------------------------
   Addieren zweier Vectoren. 
--------------------------------------------------------------------------*/
void 
add_vec (int *vec, matrix_TYP *matvec)
{
int i;
for (i=0;i<matvec->rows;i++)
   {
   vec[i] = vec[i] + matvec->array.SZ[i][0];
   }
}

/*-------------------------------------------------------------------------
    Erzeugen eines neuen Bahnelements, indem auf eine bereits bestehenden
    Wand ein Generator angewandt wird. Das ergbenis wird zurueckgegeben.  
--------------------------------------------------------------------------*/
wall_Typ *
wall_operation (derived_TYP *groupele, int groupele_ind, wall_Typ *wall, matrix_TYP *const_vec, int *genvector)
{
wall_Typ* newwall;
int* hilf_vec;
newwall=(wall_Typ*) malloc (sizeof(wall_Typ));
/* reserviert einen Speicherplatz mehr, fuer das neue Element */
newwall->product=copy_product_r(wall->product, wall->nproduct);
newwall->nproduct=wall->nproduct+1;
newwall->product[newwall->nproduct-1]=groupele_ind;
hilf_vec = subtract_vec (wall->hplane,const_vec);
newwall->hplane = multiply (groupele->element, hilf_vec);
/*printf(" dist: %d\n",newwall->dist);*/
add_vec(newwall->hplane, const_vec);
newwall->dist = distance(newwall->hplane, genvector, const_vec->rows);
free(hilf_vec);
return(newwall);
}

/*-------------------------------------------------------------------------
   Multiplikation eines Vektors mit einer Matrix. 
--------------------------------------------------------------------------*/
matrix_TYP *
mat_vec_mul (matrix_TYP *mat, int *vec)
{
matrix_TYP* hilf;
int i, j, sum;
hilf = init_mat(mat->rows,1,"k");
for (i=0;i<mat->rows;i++)
    {
    sum = 0;
    for (j=0;j<mat->cols;j++)
        {
        sum+= mat->array.SZ[i][j]*vec[j]; 
        }
    hilf->array.SZ[i][0] = sum;
    }
return(hilf);
}

/*---------------------------------------------------------------------------
   Bestimmen des Skalarproduktes zweier Vektoren (Grammatrix wird hier nicht
   benutzt.
----------------------------------------------------------------------------*/
int 
distance (int *wallvec, int *vector, int dim)
{
int sum=0;
int i,j;
/*
printf(" \n Wall:  ");
for (i=0;i<dim;i++)
printf(" %d", wallvec[i]);
printf("\n");
printf(" \n genvec:  ");
for (i=0;i<dim;i++)
printf(" %d", vector[i]);
printf("\n");
*/
for (i=0;i<dim;i++)
	{
	sum+= (wallvec[i]*vector[i]);
	}
/*
printf(" dist : %d\n", abs(sum));
*/
return(abs(sum));
}

/*-------------------------------------------------------------------------
   Prozedur zur Brechnung der Bahn unter Gruppenoperation des Vektors
   genvector. Abgespeichert wird hierbei die Zusammensetzung der
   Gruppenelemente (die auf genvector angewendet werden) als Produkt
   der Erzeuger, sowie die Gleichung die den Halbraum definiert.
   Parameter: groupgen -- Liste der Gruppenerzeuger
              genvector -- Vektor dessen Bahn berechnet werden soll
              walls -- Ausgabeparameter; Enthaelt am Ende die Berechneten
                       Halbraeume.
              grammatrix -- Eine Grammatrix des Zugrundeliegenden Raums.
-------------------------------------------------------------------------*/
void 
calc_walls (derivedsg_TYP *groupgen, int *genvector, wall_list_TYP *walls, matrix_TYP *grammatrix)
{
int start, ende;
int i, j;
int index;
matrix_TYP* const_vec;
wall_Typ* newwall;
const_vec = mat_vec_mul(grammatrix,genvector);
newwall = (wall_Typ*) malloc (sizeof(wall_Typ));
newwall->nproduct = 0;
newwall->hplane = (int*) malloc (sizeof(int)*walls->dimension);
newwall->dist = 0;
for (i=0;i<walls->dimension;i++)
    newwall->hplane[i] = 0;
insert_wall( walls, newwall );
start = 0;
ende=walls->firstfree;

do
   {
   /* deleted the next 3 lines because they don't give a clever
      orbit algorithm
   for (i=0;i<groupgen->firstfree;i++)
       {
       for (j=0;j<ende;j++) */
   for (j=0;j<walls->firstfree;j++)
       {
          /* don't bore the user */
          if (j% 100 == 0){
             fprintf(stderr," . ");
             fflush(stderr);
          }

          for (i=0;i<groupgen->firstfree;i++)
          {
          newwall=wall_operation(groupgen->list[i], i, walls->wall_list[j], const_vec, genvector);
          if (!is_in_wall_list(walls, newwall, &index))
             {
             insert_wall(walls, newwall);
             }
          else
             {
             free_wall_l(newwall);
             }
          }
       }

   /* done this way do leave joerg's frame of programming */
   start=ende=walls->firstfree;
   fprintf(stderr,"\n");
   }
while (start!=ende);
free_mat(const_vec);
}

/*-------------------------------------------------------------------------
   Ausgeben der Wanddaten. 
--------------------------------------------------------------------------*/
void 
dump_walls (wall_list_TYP *walls)
{
int i, j;
printf(" Hplanes: \n");
for (i=0; i<walls->firstfree; i++)
    {
    for (j=0; j<walls->dimension;j++)
        printf(" %d ",walls->wall_list[i]->hplane[j]);
    printf(" \n");
    }
}

/*-------------------------------------------------------------------------
   Einlesen der Grammatrix. 
--------------------------------------------------------------------------*/
matrix_TYP *
get_grammatrix (char *filename, int dim)
{
matrix_TYP* hilf;
FILE* infile;
infile = fopen(filename, "r");
hilf = fget_mat(infile);
/*dumpmat(hilf);*/
fclose(infile);
return(hilf);
}

/*
*/
void 
testrep (derivedsg_TYP *groupgen, derivedsg_TYP *generators, int dim)
{
int i,j;
matrix_TYP* hilfm;
derived_TYP* hilf;
hilf = (derived_TYP*) malloc (sizeof(derived_TYP));
hilf->nproduct = 0;
for (i=0;i<groupgen->firstfree;i++)
    {
    hilf->element = init_mat (dim, dim, "");
    for (j=0;j<dim;j++)
        hilf->element->array.SZ[j][j]=1;
    for (j=0;j<groupgen->list[i]->nproduct;j++)
        {
        hilfm = hilf->element;
        hilf->element = mat_mul(hilfm, generators->
                list[groupgen->list[i]->product[j]]->element);
        free_mat(hilfm);
        }
    if (!is_equal(groupgen->list[i], hilf))
        printf(" Representation of groupgenerator %d is wrong !\n",i);
    }
free(hilf);
}

/*
*/
void 
testrep_rel (presentation_TYP *presentation, int dim)
{
int i,j;
matrix_TYP* hilfm;
derived_TYP* hilf;
derived_TYP* eins;
hilf = (derived_TYP*) malloc (sizeof(derived_TYP));
hilf->nproduct = 0;
eins = (derived_TYP*) malloc (sizeof(derived_TYP));
eins->nproduct =0;
eins->element = init_mat (dim, dim, "");
for (j=0;j<dim;j++)
    eins->element->array.SZ[j][j]=1;
for (i=0;i<presentation->norelators;i++)
    {
    hilf->element = init_mat (dim, dim, "");
    for (j=0;j<dim;j++)
        hilf->element->array.SZ[j][j]=1;
    for (j=0;j<presentation->relators[i].lhsnproduct;j++)
        {
        hilfm = hilf->element;
        hilf->element = mat_mul(hilfm, presentation->generators->
                list[presentation->relators[i].lhsproduct[j]]->element);
        free_mat(hilfm);
        }
    if (!is_equal(eins, hilf))
        {
        printf("\n Error in relator number %d !\n",i);
        }
    }
free(hilf);
}

/*--------------------------------------------------------------------------
    Invertieren der Fundamentalbereichs-Erzeuger, da sonst die Darstellung
    als Produkt der urspruenglichen Erzeuger nicht stimmt.
--------------------------------------------------------------------------*/
void 
inverse_gen (presentation_TYP *presentation, int *orders)
{
int i;
int j;
int* ihilf;
matrix_TYP* hilf;
for (i=0;i<presentation->generators->firstfree;i++)
    {
    hilf = presentation->generators->list[i]->element;
    presentation->generators->list[i]->element = mat_inv(hilf);
    free_mat(hilf);
    ihilf = copy_product_r(presentation->generators->list[i]->product, presentation->generators->list[i]->nproduct);
    for (j=0;j<presentation->generators->list[i]->nproduct;j++)
        {
        presentation->generators->list[i]->product[j] = ihilf[presentation->generators->list[i]->nproduct-1-j];
        }
    free(ihilf);
    }
for (j=0;j<presentation->generators->firstfree;j++)
    {
    substitute_inverse(presentation, j, orders[j]); 
    }
/* kuerzen der Relatoren durch beruecksichtigen der Ordnung der Elemente */
/*
dumprelators(presentation->relators, presentation->norelators);
*/
inverse_prods(presentation, orders);
/*printf("\n after.....\n\n\n");*/
/*
dumprelators(presentation->relators, presentation->norelators);
*/
}

/*-------------------------------------------------------------------------
   Hauptprogramm. 
--------------------------------------------------------------------------*/
int 
main (int argc, char **argv)
{
int *genvector,
    *old_genvector,
    *P,
    *gen,
     tries;
matrix_TYP *gen_vector,
           *tmp,
           *graminv;
int dimension;
int order;
int found;
FILE* outdat;
FILE* infile;
wall_list_TYP* walls;
derivedsg_TYP* groupgen;
derivedsg_TYP* groupgennew;
bravais_TYP* bravgroup;
matrix_TYP* grammatrix;
fund_domain* f_domain;
presentation_TYP* presentation;
int i, j;
char filename_out[100];
char filename_in[100];
int* orders;
found = 0;

/* inserted to konform with the carat library */
read_header(argc,argv);

if (is_option('d')){
   SFLAG = 1;
}

fprintf(stderr,"\n\n Program calculating a presentation of a group\n");
fprintf(stderr," using properties of a fundamental domain.\n");
fprintf(stderr,"\n (c) 1995 Lehrstuhl B Mathematik, RWTH Aachen\n\n");

if ((is_option('h') && optionnumber('h') == 0) || FILEANZ < 1)
{
     printf("Usage:\n");
     printf("\n");
     printf(" Roundcor file1 [file2] [-f] [-g] [-C] [-G] [-M] [-d] [-h]\n");
     printf("\n");
     printf("Where file1 contains a bravais_TYP.\n");
     printf("\n");
     printf("The program calculates a presentation for the FINITE group G\n");
     printf("in file.\n");
     printf("The presentation is written to file2, and if this isn't given\n");
     printf("to stdout (all other output is written to stderr).\n");
     printf("\n");
     printf("The options are:\n");
     printf("\n");
     printf(" -f   : read a positive definite symetric G-invariant form\n");
     printf("        from a file called <file>.gram .\n");
     printf(" -g   : read a generic vector for G from <file>.gen_vector.\n");
     printf("        If this option is not present the program will try\n");
     printf("        one and exit if this one isn't generic.\n");
     printf(" -C   : Output the presentation in Cayley format on a file\n");
     printf("        called <file>.cayley.short .\n");
     printf(" -G   : Output the presentation in Gap format on a file\n");
     printf("        called <file>.gap.short .\n");
     printf(" -M   : Output the presentation in Magma format on a file\n");
     printf("        called <file>.magma.short .\n");
     printf(" -d   : DEBUGGING MODE; DO NOT USE.\n");
     printf(" -h   : Gives you this help.\n");
     printf("\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
     }

fprintf(stderr,"\n The name of the inputfile: %s\n",FILENAMES[0]);

if (!exist(FILENAMES[0]))
   {
   printf("\n Error: Filename '%s' does not exist !\n", FILENAMES[0]);
   exit(3);
   }

bravgroup = get_bravais(FILENAMES[0]);
dimension = bravgroup->dim;
order = bravgroup->order;
fprintf(stderr, "\n Order of the group: %d\n", order);
groupgen = init_derivedgen(dimension);
convert_bravais(bravgroup, groupgen);
walls = init_wall_list(dimension);

/* read the gram matrix */
if (is_option('f')){
   sprintf(filename_in, "%s.gram", FILENAMES[0]);
   grammatrix = get_mat(filename_in);
}
else{
   tmp = init_mat(bravgroup->dim,bravgroup->dim,"1");
   grammatrix = rform(bravgroup->gen,bravgroup->gen_no,tmp,101);
   free_mat(tmp);
}

free_bravais(bravgroup);

if (grammatrix->rows!=dimension)
   {
   printf("\n Error: grammatrix has wrong dimension! \n");
   exit(3);
   }

/* read the generic vector if nessesary */
tries = 0;
do{
if (!is_option('g'))
	{
        tries++;
        if (tries>30){
           fprintf(stderr,"I tried 30 vectors, all of them weren't generic.\n");
           fprintf(stderr,"Try again or feed me one.\n");
           exit(3);
        }
        else if(tries>1){
           fprintf(stderr,"Trying the %d-th random vector\n",tries);
           free_wall_list(walls);
           free(genvector);
           free(old_genvector);
           walls = init_wall_list(dimension);
        }

	genvector = init_vector(dimension);
	old_genvector = init_vector(dimension);
	set_vector(old_genvector, dimension,tries);

        /* changed to conform to the way the read vector is treated,
           (ie. multiply it with the inverse of the gram matrix
           tilman 25/3/97: */
        graminv = mat_inv(grammatrix);
        for (i=0;i<dimension;i++)
           for (j=0;j<dimension;j++)
              genvector[i] += graminv->array.SZ[i][j] * old_genvector[j];

        free_mat(graminv);

        if (is_option('d')){
           printf("%dx1 %% genvector\n",dimension);
           for (i=0;i<dimension;i++) printf("%d \n",genvector[i]);
           printf("%dx1 %% old_genvector\n",dimension);
           for (i=0;i<dimension;i++) printf("%d \n",old_genvector[i]);
        }
	}
else
	{
	sprintf(filename_in, "%s.gen_vector", FILENAMES[0]);
	infile = fopen(filename_in,"r");

        /* changed the next line to the next 13 25/3/97 tilman, altough
           it might look stupid, it is the easiest way to implement it
           and leave the rest of the program unchanged
	gen_vector = fget_mat(infile); */

	tmp = fget_mat(infile);
        /* check trivialities */
        if (tmp->rows != dimension || tmp->cols !=1){
           fprintf(stderr,"The generic vector should be %dx1\n",dimension);
           exit(3);
        }
        graminv = mat_inv(grammatrix);
        gen_vector = mat_mul(graminv,tmp);
        free_mat(graminv);
        rat2kgv(gen_vector);
        gen_vector->kgv = 1;
        Check_mat(gen_vector);

        /* we will need the not transformed generic vector again */
	old_genvector = init_vector(dimension);
	for (i=0;i<dimension;i++)
		old_genvector[i] = tmp->array.SZ[i][0];
        free_mat(tmp);

	fclose(infile);

        if (is_option('d')){
           printf("\n");
	   printf("Read generic vector: (in this case gram^(-1)*genvector)\n");
	   put_mat(gen_vector,NULL,NULL,2);
        }

	genvector = init_vector(dimension);
	for (i=0;i<dimension;i++)
		{
		genvector[i] = gen_vector->array.SZ[i][0];
		}
	free_mat(gen_vector);
	}

/* Generieren der "Waende" (Bahnenalgorithmus) */
fprintf(stderr," Generating the walls....\n");
calc_walls(groupgen, genvector, walls, grammatrix);

} while(walls->firstfree<order && !is_option('g'));

if (walls->firstfree<order){
   fprintf(stderr," Error in roundcor: stabelizer of the generating vector\n");
   fprintf(stderr," is bigger than one!\n");
   fprintf(stderr," Calculated Walls: %d\n Order of the group: %d\n",
           walls->firstfree, order);
   exit(3);
}

/*dump_walls ( walls );*/

fprintf(stderr," Order of the group: %d \n", walls->firstfree);
/* Bestimmung des Fundamentalbereichs */
f_domain = fundamental_domain( walls, groupgen );
free_wall_list(walls);
/***************
for (i=0;i<f_domain->wall_no;i++)
   {
   printf(" \n\n Wall associate matrix: \n");
   dumpmat(f_domain->wall[i]->mat);
   printf("\n Wallequation: ");
   for (j=0;j<dimension;j++)
      printf(" %d ",f_domain->wall[i]->gl[j]);
   printf("\n");
   }
***************/
/**************
printf(" Vertices: \n");
for (i=0;i<f_domain->vert_no;i++)
    {
    for (j=0;j<dimension;j++)
        printf(" %d ",f_domain->vert[i]->v[j]);
    printf("\n");
    }
**************/

/* changed 25/3/97 tilman because it causes trouble in free, and I didn't
find a better solution right now
orders = (int*) malloc (sizeof(int)*f_domain->wall_no); */
orders = (int*) calloc (f_domain->wall_no+groupgen->firstfree,sizeof(int));

presentation = generate_pres (f_domain, orders);

/* output for debugging */
if (is_option('d')){
   sprintf(filename_out, "%s.gap.1", FILENAMES[0]);
   printf(" 2. Input for GAP : %s\n", filename_out);
   gap_out(presentation, filename_out);
   if (optionnumber('d') == 1) exit(1);
}


/*
printf("\n Testing the relators after generate_pres... \n\n");
dumprelators(presentation->relators, presentation->norelators);
testrep_rel (presentation, dimension);
sprintf(filename_out, "%s.cayley.test0", FILENAMES[0]);
cayley_out(presentation, filename_out);
*/
/*
cayley_out(presentation);
exit();
*/
inverse_gen(presentation, orders);
/*
printf(" testing presentation->generators ...\n");
testrep (presentation->generators, groupgen, dimension);
*/
groupgennew = copy_derivedsg (presentation->generators);
old_as_new(groupgen, groupgennew, dimension);
/*
printf("\n Testing groupgen ...\n");
testrep (groupgen, groupgennew, dimension);
printf("\n Testing groupgennew ...\n");
testrep (groupgennew, groupgen, dimension);
*/
/*
for (i=0;i<groupgen->firstfree;i++)
    {
    printf(" Oldgen no %d : \n",i);
    dumpmat(groupgen->list[i]->element);
    dumpproduct(groupgen->list[i]->product,groupgen->list[i]->nproduct);
    }
for (i=0;i<groupgennew->firstfree;i++)
    {
    printf(" newgen no %d : \n",i);
    dumpmat(groupgennew->list[i]->element);
    dumpproduct(groupgennew->list[i]->product,groupgennew->list[i]->nproduct);
    }
dumplist(presentation->generators);
*/
tietzetrans(groupgen, groupgennew, presentation);
/*
printf("\n Testing the relators after Tietzetrans... \n\n");
testrep_rel (presentation, dimension);
*/
if (is_option('d')) fprintf(stderr,"\n Output files are: \n");

if (is_option('C') && is_option('d')){
   sprintf(filename_out, "%s.cayley", FILENAMES[0]);
   printf(" 1. Input for Cayley-Software : %s\n", filename_out);
   cayley_out(presentation, filename_out);
}

if (is_option('G') && is_option('d')){
   sprintf(filename_out, "%s.gap", FILENAMES[0]);
   printf(" 2. Input for GAP : %s\n", filename_out);
   gap_out(presentation, filename_out);
}

if (is_option('d')){
   sprintf(filename_out, "%s.carat", FILENAMES[0]);
   printf(" 3. Input for CARAT : %s\n", filename_out);
   carat_out(presentation, filename_out);
}

/* this program is not supported anymore */
if (FALSE){
   sprintf(filename_out, "%s.rel", FILENAMES[0]);
   printf(" 4. Input for ZASS : %s\n", filename_out);
   outdat = fopen(filename_out,"w");
   dumprelators_file(outdat, presentation->relators, presentation->norelators);
}

/* inserted for debugging */
if (is_option('d')){
   sprintf(filename_out, "%s.gap.1", FILENAMES[0]);
   printf(" 2. Input for GAP : %s\n", filename_out);
   gap_out(presentation, filename_out);
   if (optionnumber('d') == 1) exit(1);
}

/* unterrelatoren werden gesucht und eliminiert */
fprintf(stderr," Splitting relators into subrelators ... \n");
fprintf(stderr," no of relators before: %d \n",presentation->norelators);

sub_relator_search(presentation,old_genvector, dimension);

/* inserted for debugging */
if (is_option('d')){
   sprintf(filename_out, "%s.gap.2", FILENAMES[0]);
   printf(" 2. Input for GAP : %s\n", filename_out);
   gap_out(presentation, filename_out);
   if (optionnumber('d') == 2) exit(2);
}

/* Duplikate werden eliminiert */
calc_orders(presentation, orders, dimension);
fprintf(stderr," Eliminating duplicate relators ... \n");
shorten_presentation(presentation, orders);
free(orders);

/* changed tilman 26/3/97, reorganized the output */
fprintf(stderr,"\n Output files are: \n");

if (is_option('C')){
   sprintf(filename_out, "%s.cayley.short", FILENAMES[0]);
   fprintf(stderr," Input for Cayley-Software : %s\n", filename_out);
   cayley_out(presentation, filename_out);
}

if (is_option('G')){
   sprintf(filename_out, "%s.gap.short", FILENAMES[0]);
   fprintf(stderr," Input for GAP : %s\n", filename_out);
   gap_out(presentation, filename_out);
}

/* reorganized output tilman 15/5/97
sprintf(filename_out, "%s.carat.short", FILENAMES[0]);
fprintf(stderr," 3. Input for CARAT : %s\n", filename_out); */
if (FILEANZ > 1)
   carat_out(presentation, FILENAMES[1]);
else
   carat_out(presentation, NULL);

if (is_option('M')){
   sprintf(filename_out, "%s.magma.short", FILENAMES[0]);
   fprintf(stderr," 4. Input for MAGMA : %s\n", filename_out);
   magma_out(presentation, filename_out);
}

/* the program is not supported anymore */
if (FALSE){
   sprintf(filename_out, "%s.rel.short", FILENAMES[0]);
   printf(" 5. Input for ZASS : %s\n", filename_out);
   outdat = fopen(filename_out,"w");
   dumprelators_file(outdat, presentation->relators, presentation->norelators);
   fclose(outdat);
}

/*
free_derivedsgcomp(groupgennew);
free_presentation(presentation); 
free_fund_domain (f_domain);
free_vector(genvector);
free_derivedsgcomp(groupgen);
free_mat(grammatrix);
*/

if (is_option('d')){
   fprintf(stderr," SCOUNT: %d\n",SCOUNT);
   pointer_statistics(0,0);
}


exit(0);
}
