#include "typedef.h"
#include "getput.h"
#include "dertype.h"
#include "tietzetrans.h"

extern int* copy_product_r (int*, int);
extern matrix_TYP* mat_inv (matrix_TYP*);
extern matrix_TYP* copy_mat (matrix_TYP*);
extern matrix_TYP* mat_mul (matrix_TYP*, matrix_TYP*);
extern derivedsg_TYP* init_derivedgen (int);
extern void free_mat(matrix_TYP*);
extern derived_TYP* copy_derived(derived_TYP*);

/* lokale Typdefinitionen */
typedef struct{
            int w[2];
            int* v;
            int nv;
            } vert_list_TYP;

/* Prototypenliste */
int first_free_ar( int*, int);
void insert_ar(int*, int, int);
vert_list_TYP* generate_vert_list (fund_domain*, int, int*);
int operation_on_vert (int, vert_list_TYP*, int, fund_domain*, int,int);

/*------------------------------------------------------------------------------
 Funktion, die in einem Array die erste Null sucht und Stelle zurueckgibt 
 Parameter:
     ar - Feld von Integer (Eingabe)
     ar_len - Laenge des Feldes (Eingabe)
 Rueckgabewert:
     Position der ersten Null im Feld.
   oder
     -1, falls keine Null vorhanden ist.
-----------------------------------------------------------------------------*/
int first_free_ar (ar, ar_len)
int* ar;
int ar_len;
{
int i;
for (i=0;i<ar_len;i++)
   if (!ar[i]) return (i);
return(-1);
}

/*----------------------------------------------------------------------------
   Funktion, die an der Stelle no ueberprueft, ob eine Eins gesetzt ist und
   ggf. eine setzt. 
  Parameter:
     ar - Integerfeld (Ein/Ausgabe)
     no - Position im Feld, die gesetzt werden soll. (Eingabe)
     ar_len - Laenge des Integerfeldes ar. (Eingabe)
----------------------------------------------------------------------------*/
void insert_ar (ar, no, ar_len)
int* ar;
int ar_len;
int no;
{
if ((no>-1)&&(no<ar_len) && (!ar[no]))
   ar[no]=1;
}

/*--------------------------------------------------------------------------
    Gruppenelement g auf den Integer-Vektor mit dem Index
    vert in der Liste fudo->vert[] anwenden.
    Parameter:
              g 	- Matrix des Gruppenelements
              vert 	- Index des Vektors in der Liste fudo->vert[]
              fudo 	- Fundamentalbereich, enthaelt unter anderem den
                   	  Vektor, auf den die Operation angewendet werden soll.
              end 	- Liste mit den Zielecken.
              endlen 	- Anzahl der Zielecken.
              dim 	- Dimension der Raumes.

    Rueckgabewert:
    Die Funktion liefert 1 (true) wenn die Ecke durch g auf eine der Ecken
    in dem Feld end abgebildet wird, sonst Null.  
    
----------------------------------------------------------------------------*/
int op_vector_normal ( g, vert, fudo, end, endlen, dim)
matrix_TYP* g;
int vert;
fund_domain* fudo;
int* end;
int endlen;
int dim;
{
int erg;
int i,j;
int is_eq;
int* v;
int k;
int is_in_list;
v=(int*) malloc (sizeof(int)*dim);
for (i=0;i<dim;i++)
    {
     v[i]=0;
     for (j=0;j<dim;j++)
         {
         v[i]+= (g->array.SZ[j][i]*fudo->vert[vert]->v[j]);
         }
    }
is_in_list = 0;
for (k=0;(k<endlen)&&(!is_in_list);k++)
   {
   is_eq = 1;
   for (j=0;j<dim;j++)
       is_eq = is_eq && (v[j]==fudo->vert[end[k]]->v[j]);
   is_in_list = is_eq || is_in_list;
   }
free(v);
return(is_in_list);
}

/*--------------------------------------------------------------------------
    Inverses des Gruppenelements g auf den Integer-Vektor mit dem Index
    vert in der Liste fudo->vert[] anwenden.
    Parameter:
              g - Matrix des Gruppenelements
              vert - Index des Vektors in der Liste fudo->vert[]
              fudo - Fundamentalbereich, enthaelt unter anderem den
                     Vektor, auf den die Operation angewendet werden soll.
              dim - Dimension der Raumes.
   Rueckgabewert:
   Nummer der vert, auf die vert abgebildet wird. 
   -1, wenn kein vert identifiziert.
----------------------------------------------------------------------------*/
int op_vector ( g, vert, fudo, dim)
matrix_TYP* g;
int vert;
fund_domain* fudo;
int dim;
{
int erg;
int i,j;
int is_eq;
matrix_TYP* hilf;
int* v;
v=(int*) malloc (sizeof(int)*dim);

/* hilf = mat_inv( g ); */
hilf = copy_mat( g );

/* output for debugging */
if (is_option('d')){
   printf(" hilf: \n");
   put_mat(hilf,NULL,NULL,2);
}

for (i=0;i<dim;i++)
    {
     v[i]=0;
     for (j=0;j<dim;j++)
         {
         v[i]+= (hilf->array.SZ[j][i]*fudo->vert[vert]->v[j]);
         }
    }
/* suche nun den Index der berechneten Vertex */
for (i=0;i<fudo->vert_no;i++)
    {
    is_eq=1;
    for (j=0;j<dim;j++)
        is_eq = is_eq && (v[j]==fudo->vert[i]->v[j]);
    if (is_eq)
        {
        free(v);
        free_mat(hilf);
        return(i);
        }
    }
free(v);
free_mat(hilf);
return(-1);
}

/*-----------------------------------------------------------------------
   Operation auf einer Wand ausfuehren (d.h. auf allen verts, die auf 
   dieser Wand liegen.
   Parameter:
          wall - Wand, um die es geht.
          fudo - Fundamentalbereich mit allen Daten ueber verts und Waende
          dim - Dimension des Raumes.
   Rueckgabewert:
          Nummer der resultierenden Wand oder
          -1, wenn die Ergebniswand nicht identifiziert werden konnte.
------------------------------------------------------------------------*/
int operation_on_wall( wall, fudo, dim, max_wallno)
int wall;
fund_domain* fudo;
int dim;
int max_wallno;
{
int i, k, j, l;
int* v;
int erg;
int is_in_list, is_wall_j;
v = (int*) malloc (sizeof(int)*10000);
/*
v = (int*) malloc (sizeof(int)*max_wallno);
v = (int*) malloc (sizeof(int)*fudo->vert[0]->wall_no);
*/
/* Auf alle vertices, die auf der angegebenen Wand liegen (angegeben ist
   der Index der Wand in der Liste wall, welches ein Attribut des 
   Fundamentalbereichs ist), wende die Matrix an, die mit der Wand wall
   assoziiert ist.
*/
k=0;
/* ueber alle verts */
for (i=0; i<fudo->vert_no;i++)
    {
    is_in_list=0;
   /* ueber alle Waende auf denen dieser vert liegt */
    for (j=0;(j<fudo->vert[i]->wall_no) && (!is_in_list);j++)
        {
        is_in_list = is_in_list || (wall==fudo->vert[i]->wall[j]);
        }
  /* vert liegt auf der Wand wall */
    if (is_in_list)
       {
     /* inverses der mit wall assoziierten Matrix auf den zur vert i
        gehoerenden Vector anwenden. Nummer der Ergebnisvert abspeichern */
       v[k] = op_vector( fudo->wall[wall]->mat, i, fudo, dim );
  /*printf(" wall: %d, vertice: %d, ergvert: %d\n", wall, i , v[k]);*/
       if (v[k]<0) 
          printf("\n Error in op_vector: result is not a vert of the fundamental domain !\n\n");
       k++;
       } 
    }
/* identifiziere nun aus der Liste der vertices v[] um welche Wand
   es sich handelt. Dabei untersuche alle Vertices in v[] und schaue,
   welche Wand auf allen draufliegt (muss eindeutig sein).
*/
for (j=0;j<fudo->wall_no;j++)
    {
    is_wall_j=1;
    for (i=0;i<k;i++)
       {
       is_in_list = 0;
       for (l=0;(l<fudo->vert[v[i]]->wall_no)&&(!is_in_list);l++)
           is_in_list = is_in_list || (fudo->vert[v[i]]->wall[l]==j); 
       is_wall_j= is_wall_j && is_in_list;
       }
    if (is_wall_j)
       {
       free(v);
       return(j);
       }
    }
/* keine Wand identifiziert: */
free(v);
return(-1);
}

/*---------------------------------------------------------------------------
   Operation auf einer Ecke Durchfuehren.   
   Dabei wird auf die Wand, die ueberquert werden soll, das Inverse
   des assoziierten Gruppenelements angewendet (operation_on_wall). 
   Das Ergebnis ist eine andere Wand im Fundamentalbereich. 
   Um die Ecke zu ermitteln
   betrachtet man nun alle Ecken, die auf dieser Wand liegen.
   Die Ecke, die durch Operation mit dem Gruppenelement auf die
   Startecke abgebildet wird, ist die gesuchte.
   Das Ergebnis der Operation ist eine Ecke; deren Index in der 
   vert_list wird als Ergbenis zurueckgegeben.
   Parameter:
       act_vert - Index der aktuell zu behandelnden Ecke in der vert_list
       vert_list - Liste aller Ecken des Fundamentalbereichs
       vert_no - Anzahl Ecken in vert_list
       fudo - Der Fundamentalbereich
       dim - Dimension des Raumes in dem gerechnet wird.
   Rueckgabewert:
       resultierende Ecke
---------------------------------------------------------------------------*/
int operation_on_vert (act_vert, vert_list, vert_no, fudo, dim, maxwalls)
int act_vert;
vert_list_TYP* vert_list;
int vert_no;
fund_domain* fudo;
int dim;
int maxwalls;
{
int wall_1;
int erg;
int i;
int j;
int l;
int k;
int* erg_list;
int erg_no;
int is_eq;
/* Die neue Wand bestimmen, die unter der Operation, des Inversen des mit
   der Wand 0 assoziierten Gruppenelements, entsteht. 
   Wand 0, weil:
   Beim Start ist es egal welche Wand genommen wird und in den nachfolgenden
   Iterationen wird falls noetig vorher eine Vertauschung vorgenommen.
*/
wall_1 = operation_on_wall (vert_list[act_vert].w[0],fudo, dim, maxwalls);
if (wall_1<0) 
   printf("\n Error in operation_on_wall: Couldn't identify wall !!\n");
/* bestimme welche Ecken durch wall_1 mit anderen Waenden gebildet werden,
   Moeglichkeiten werden in erglist abgespeichert. */
k=0;
/*
erg_list = (int*) malloc (sizeof(int)*dim);
*/
erg_list = (int*) malloc (sizeof(int)*10000);
for (i=0;i<vert_no;i++)
    {
    if ((vert_list[i].w[0]==wall_1) || (vert_list[i].w[1]==wall_1))
       {
       erg_list[k]=i;
       k++;
       }
    }
erg_no = k;
/* Stelle fest welche der Ecken auf die Startecke abgebildet wird
   unter der Operation des Gruppenelements das mit der Wand 0 identifiziert
   wird.
*/
/* stelle nun fest welche der Ecken in erg_list auf die Startecke act_vert
   abgebildet wird */
/* Dazu werden zunaechst die verts bestimmt, die auf den Ecken in erg_list
   liegen. */
erg=0;
for (i=0;(i<erg_no)&&(erg==0);i++)
    {
/* Operation auf alle verts der Ecke anwenden (jede vert muss auf eine
   vert der Ecke act_vert (abgespeichert in on_walls) gehen, wenn es die
   gesuchte ist. */
    is_eq = 1;
    for (j=0;(j<vert_list[erg_list[i]].nv)&&(is_eq);j++)
        {
        is_eq = op_vector_normal(fudo->wall[vert_list[act_vert].w[0]]->mat,
                vert_list[erg_list[i]].v[j], fudo, vert_list[act_vert].v,
                vert_list[act_vert].nv, dim); 
        }
    if (is_eq) 
       {
       erg = i;
       }
    }
erg = erg_list[erg];
/* Falls noetig Waende vertauschen. */
if (vert_list[erg].w[0]==wall_1)
   {
   vert_list[erg].w[0]=vert_list[erg].w[1];
   vert_list[erg].w[1]=wall_1; 
   }
free(erg_list);
return(erg);
}

/*----------------------------------------------------------------------------
  Speicherplatz der Ecken-liste wieder freigeben. Wird in generate_vert_list
  angefordert.
----------------------------------------------------------------------------*/
void free_vert_list( list, length) 
vert_list_TYP* list;
int length;
{
int i;
for (i=0;i<length;i++)
    {
    free(list[i].v);
    }
free(list);
}

/*------------------------------------------------------------------------------
    Erstellen einer Liste mit allen Ecken (i.e. Schnitte von Waenden, i.e.
    Raeumen der codim 2).Dabei sind Waende als Raeume der codim 1 zu verstehen.
    Parameter:
            fudo - Fundamentalbereich
            dim - Dimension des Raumes
            no - Rueckgabewert; Anzahl der generierten Ecken
    Rueckgabewert:
            Liste der generierten Ecken.
-----------------------------------------------------------------------------*/
vert_list_TYP* generate_vert_list (fudo, dim, no)
fund_domain* fudo;
int dim;
int* no;
{
int i,j,k,l,m;
int jok, iok;
int count;
int vert_anz;
vert_list_TYP* erg;
vert_anz = (fudo->wall_no-1)*fudo->wall_no/2;
erg = (vert_list_TYP*) malloc (sizeof(vert_list_TYP)*vert_anz);
count = 0;
m=0;
for (i=1; i<fudo->wall_no;i++)
    {
    for (j=0;j<i;j++)
        {
        erg[count].w[0]=j;
        erg[count].w[1]=i;
        erg[count].v = (int*) malloc (sizeof(int)*10000);
        for (k=0;k<fudo->vert_no;k++)
            {
            jok=0;
            iok=0;
            for (l=0;l<fudo->vert[k]->wall_no;l++)
                {
                jok = jok || (fudo->vert[k]->wall[l]==j);
                iok = iok || (fudo->vert[k]->wall[l]==i);
          /*      if (fudo->vert[k]->wall[l]==j)
           **         jok=1;*/
          /*    if (fudo->vert[k]->wall[l]==i)
           **         iok=1;*/
                }                   
            if (iok && jok)
               {
               erg[count].v[m] = k;
               m++;
               }
            }
            erg[count].nv=m;
/* Ueberpruefen, ob der Schnitt der Waende einen codim-2 Raum bildet.
   (als Fundamentalbereichsbegrenzung, sonst sicherlich abgesehen von
    Parallelitaet gegeben)
*/
            if (m>=(dim-2))
               {
               /* we allocated far too much memory, so inserted this
                  22/4/97 */
               erg[count].v = (int *) realloc(erg[count].v,m * sizeof(int));

               count++;
               m=0;
               } 
            else
               {
               m=0;
               }
        }
    }
*no = count;
return (erg);
}



/*--------------------------------------------------------------------------
   Generieren einer Praesentation anhand eines Fundamentalbereichs
   Parameter:
            fudo - Fundamentalbereich, mit allen Daten ueber die 
                   begrenzenden Waende, verts etc.
   Rueckgabewert:
            Praesentation, die allerdings noch zwei Nachteile hat (mind.)
            1. Es fehlen die Generatorordnungen
            2. Die Praesentation liegt in den "falschen Generatoren" vor
 (in den Generatoren, die den Fundamentalbereich bestimmen).
--------------------------------------------------------------------------*/
presentation_TYP* generate_pres(fudo, orders)
fund_domain* fudo;
int *orders;
{
int dim;
int i;
int k;
int l;
int vert_path_len;
int act_vert, start_vert;
int vert_no;
int startwall;
int is_inverse;
int *inverse_list_i;
int ninverse_list_i;
int jk;
int maxwalls;
int *gen_list;
int gen_count;
vert_list_TYP* vert_list;
/* Feld mit Einsen und Nullen, Eins heisst, dass die entsprechende Ecke
   bereits durchlaufen wurde, Null bedeutet entsprechend das Gegenteil */
int* vert_path;
derived_TYP* hilf2;
int index;
int ende;
derived_TYP* hilf;
derivedsg_TYP* inverse_list;
presentation_TYP* presentation;
dim = fudo->wall[0]->dim;
maxwalls=0;
for (i=0;i<fudo->vert_no;i++)
    {
    if (fudo->vert[i]->wall_no>maxwalls)
        {
        maxwalls = fudo->vert[i]->wall_no;
        }
    }
/* Generieren einer Liste von "Ecken", die aus jeweils zwei Waenden bestehen
   und den vertices, die auf diesen Ecken liegen. Ecke ist dann als Schnitt
   von zwei Waenden zu verstehen. Es wird ausserdem geprueft, ob die Ecken
   einen Raum der Co-dim. 2 bilden, falls nicht werden Sie bereits aus-
   sortiert. Anzahl der gefundenen Ecken wird ebenfalls zurueckgegeben.
*/
vert_list = generate_vert_list (fudo, dim, &vert_no);
/* Speicherplatz anfordern, fuer die Liste der bereits durchlaufenen
   Ecken. */
vert_path = (int*) malloc (sizeof(int)*vert_no);
for (i=0;i<vert_no;i++)
    vert_path[i]=0;
/* Praesentation initialisieren */
presentation = (presentation_TYP*) malloc (sizeof(presentation_TYP));
/* Generatoren definieren, sind bei den Waenden abgespeichert
   Problem: Wenn Inverses eines Generators auch eine Wand des Fundamental-
   bereichs ist, soll diese Wand natuerlich nicht uebernommen werden. */
presentation->generators = init_derivedgen(dim);
inverse_list = init_derivedgen(dim);
inverse_list_i = (int*) malloc (sizeof(int)*fudo->wall_no);
ninverse_list_i = 0;
gen_list = (int*) malloc (sizeof(int)*fudo->wall_no);
gen_count = 0;
presentation->norelators = 0;
presentation->relators = (relator_TYP*) malloc (sizeof(relator_TYP)*EXT_SIZE);
presentation->ext_factor = 1;

for (i=0;i<fudo->wall_no;i++)
    {
    hilf = (derived_TYP*) malloc (sizeof(derived_TYP));
    hilf->nproduct =0;
    hilf->element = copy_mat (fudo->wall[i]->mat);
    hilf2 = (derived_TYP*) malloc (sizeof(derived_TYP));
    hilf2->nproduct = 0;
    hilf2->element = mat_inv (fudo->wall[i]->mat);
/* Weder das Element noch das Inverse des Elements sind bereits in der
   Generatorliste */
    if ((!is_in_list(presentation->generators, hilf2, &index))&&
      (!is_in_list(presentation->generators, hilf, &index)))
       {
       free_mat(hilf2->element);
       free(hilf2);
       hilf->nproduct = fudo->wall[i]->nproduct;
       hilf->product = copy_product_r( fudo->wall[i]->product, hilf->nproduct);
       insert_list(presentation->generators, hilf);
       calc_orders_s(presentation, orders, dim);
       gen_list[i]=gen_count;
       gen_count++;
       }
    else
       {
       is_in_list(presentation->generators, hilf2, &index);
       free_mat(hilf2->element);
/* Nummer der Wand speichern, die in der Inverselist zu finden ist */
       inverse_list_i[ninverse_list_i]=i;
       ninverse_list_i++;
/* Darstellung des Inversen als Product der Generatoren bestimmen */
       hilf->nproduct=orders[index]-1;
       hilf->product = (int*) malloc (sizeof(int)*hilf->nproduct);
       for (l=0;l<hilf->nproduct;l++)
           {
           hilf->product[l]=index;
           }
       hilf2->element = copy_mat(hilf->element);
       insert_list(inverse_list, hilf);
/* War nur zur Ueberbrueckung notwendig
       hilf2->nproduct = fudo->wall[i]->nproduct;
       hilf2->product = copy_product_r( fudo->wall[i]->product, hilf2->nproduct);
       insert_list(presentation->generators, hilf2);
*/
       }
    }
/* Wenn alle Ecken einmal durchlaufen, dann ist die Praesentation fertig
   (bis auf die Ordnung der Gruppenerzeuger, die extra brechnet werden) */
/*****
printf (" \n Waende mit inversem im Fundi \n ");
dumpproduct(inverse_list_i, ninverse_list_i);
****/
vert_path_len = 0;

fprintf(stderr," The set of generators derived from the fundamental");
fprintf(stderr," domain has cardinality: %d\n",
                    presentation->generators->firstfree);

while (vert_path_len<vert_no)
     {
/*****
     printf(" presentation->norelators %d \n",presentation->norelators);
*****/
/* Um die Ecken  "herumwandern" */
/* Suche erste Ecke, die nicht im Pfad steht, diese Ecke wird dann
   "umkreist". Startecke wird besetzt. */
     start_vert = first_free_ar(vert_path,vert_no);
     if (start_vert==(-1))
         {
         printf("\n\n Error in Generate_pres: wrong start_vert !!!\n\n");
         exit(3);
         }
/* Gehe solange weiter bis die Startecke wieder erreicht wurde */      
     act_vert = start_vert;
     startwall = vert_list[start_vert].w[0];
/* neuen Relator initialisieren */
     if (presentation->norelators>=(presentation->ext_factor*EXT_SIZE))
        {
        presentation->ext_factor++;
        presentation->relators=(relator_TYP*) realloc (presentation->relators,
            sizeof(relator_TYP)*presentation->ext_factor*EXT_SIZE);
        }
   /* Schreibabkuerzung */
     k=presentation->norelators;
     presentation->relators[k].lhsnproduct = 0;
     presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int));
     presentation->relators[k].rhsnproduct = 0;
     l=0;

     do
       {
       /* Einfuegen in den Pfad, falls die Ecke noch nicht durchlaufen wurde */
       insert_ar(vert_path, act_vert, vert_no);
       /* Gruppenelement, das mit der Wand 0 der aktuellen Ecke act_vert
          assoziert ist wird an den Relator drangehaengt. */
       presentation->relators[k].lhsnproduct++;
       presentation->relators[k].lhsproduct = (int*) realloc 
            (presentation->relators[k].lhsproduct, sizeof(int)*
               presentation->relators[k].lhsnproduct);
/* Beachte: w[0] ist eine Wandnummer, also muss eine Transformation auf
   Generatornummer erfolgen und dann das Inverse in das Produkt eingefuegt
   werden. */
       presentation->relators[k].lhsproduct[
           presentation->relators[k].lhsnproduct-1]
           = vert_list[act_vert].w[0];
/* Umsetzung der Wandnummer in Wort von Generatornummern */
       is_inverse=0;

       for (jk=0;jk<ninverse_list_i;jk++)
           {
           if (inverse_list_i[jk]==vert_list[act_vert].w[0])
              {
              insert_in_product(&presentation->relators[k].lhsproduct, presentation->relators[k].lhsnproduct-1, &presentation->relators[k].lhsnproduct,inverse_list->list[jk]->product, inverse_list->list[jk]->nproduct);
              ende = is_identity_new(presentation,inverse_list->list[jk]->nproduct, dim);
              is_inverse=1;
              }
           }
       if (!is_inverse)
          {
           presentation->relators[k].lhsproduct[
           presentation->relators[k].lhsnproduct-1]
           = gen_list[vert_list[act_vert].w[0]];
          ende = is_identity_new(presentation,1, dim);
          }
       act_vert = operation_on_vert(act_vert, vert_list, vert_no, fudo, dim,maxwalls);
       l++;
       if (l>3000) 
          {
          if (l % 20 == 0) printf("l= %d\n",l);
          printf("is_inverse %d\n",is_inverse);

          printf("Warning in generate_pres: More than 3000 Elements in a\n");
          printf("Relator. Program terminates now to avoid a possible");
          printf("infinite loop.\n");
          exit(3);
          }
       }
       while (!ende);
/*     while ((act_vert!=start_vert)||(vert_list[act_vert].w[0]!=startwall));*/
  /* is_identity prueft, ob die Relation norelators (also die letzte in der
     Liste), konsistent ist (d.h. lhsproduct ist Identitaet). */
       /*while (!is_identity(presentation, dim));*/
/* Berechnung der neuen Pfadlaenge */
     vert_path_len=0;
     for (i=0;i<vert_no;i++)
         vert_path_len+=vert_path[i];
/* */
     presentation->norelators++;
     }
free(vert_path);
free_vert_list(vert_list, vert_no);
free(inverse_list_i);
free_derivedsg(inverse_list);
return(presentation);
} 
