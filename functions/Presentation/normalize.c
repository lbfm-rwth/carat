#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"

#define MAX_GEN 24

void shift_cyclic_left(char*, int, int);
void delete_relator_s (char**, int*, int, int);
int eliminate_duplicates (char**,int*,int);
void normalize_single ( char*, int, char*);
void normalize (char**, int*,int, char*);
char* inverse_relator (char*, int, int*, char*);
int eliminate_inverse( char**, int*, int, int*, char*);
int index_table(char*,int);





/*------------------------------------------------------------------------
  Prozedur, die einen String (Character) zyklisch verschiebt und zwar um
  no_steps nach links.
  Parameter:
        string_ele - Input und Output - string, der den zu verwschiebenden
                     String enthaelt.
        string_length - Laenge des Strings string_ele
        no_steps - Anzahl der Stellen, um die Verschoben werden soll.
-------------------------------------------------------------------------*/
void shift_cyclic_left(string_ele,string_length, no_steps)
char* string_ele;
int string_length;
int no_steps;
{
char* mem;
int i;
mem = (char*)malloc (sizeof(char)*string_length);
if (no_steps>=string_length)
   {
   no_steps = no_steps % string_length;
   }
for (i=0;i<string_length;i++)
    {
    mem[i]=string_ele[(i+no_steps)%string_length];
    }
for (i=0;i<string_length;i++)
    {
    string_ele[i]=mem[i];
    }
free(mem);
}

/*----------------------------------------------------------------------------
   Prozedur, die einen Relator loescht.
   Parameter:
      relators_string - Liste der Relatoren
      relators_length - Liste der Laengen der Relatoren
      relator_count - Anzahl der Relatoren
      relatorindex - Index des zu loeschenden Relators
   Beachte: Die Anzahl der Relatoren muss separat aktualisiert werden.
----------------------------------------------------------------------------*/
void delete_relator_s (relators_string, relators_length, relator_count,
                       relatorindex)
char** relators_string;
int* relators_length;
int relator_count;
int relatorindex;
{
free(relators_string[relatorindex]);
if ((relator_count-1)==relatorindex)
   {
    relators_length[relatorindex]=0;
    relators_string[relator_count-1] = NULL;
   }
else
   {
   relators_string [relatorindex] = relators_string[relator_count-1];
   relators_string[relator_count-1] = NULL;
   relators_length[relatorindex] = relators_length[relator_count-1];
   relators_length[relator_count-1] = 0;
   }
}

/*--------------------------------------------------------------------------
   Funktion, die alle doppelt vorkommenden Relatoren eliminiert.
   Parameter:
       relators_string - Liste der Relatorenstrings
       relators_length - Liste der Laengen der Relatoren
       relator_count - Anzahl der Relatoren
   Rueckgabewert:
       Aktualisierte Anzahl der Relatoren.
-------------------------------------------------------------------------*/
int eliminate_duplicates (relators_string, relators_length, relator_count)
char** relators_string;
int* relators_length;
int relator_count;
{
int is_equal;
int i;
int j;
int k;
for (i=0;i<relator_count-1;i++)
    {
    is_equal =0;
    for (j=i;j<relator_count-1;)
        {
        if (!is_equal)
           {
           j++;
           }
        else
           {
           is_equal = 0;
           }
        if (relators_length[i]==relators_length[j])
           {
           is_equal = 1;
           for (k=0;(k<relators_length[i])&&(is_equal);k++)
               {
                if (relators_string[i][k]!=relators_string[j][k])
                   {
                    is_equal =0;
                   }
               } 
           if (is_equal)
               {
               delete_relator_s(relators_string, relators_length, relator_count,j);
               relator_count--;
               }
           }
        }
    }
return (relator_count);
}

/*----------------------------------------------------------------------------
  Funktion, die den Index eines Characters in einer Tabelle liefert.
  Parameter:
    char_table - Charactertabelle
    character - Character, der in der Tabelle gesucht werden soll.
  Rueckgabewert:
    Ermittelter Index
----------------------------------------------------------------------------*/
int index_table(char_table,character)
char* char_table; 
int character;
{
int i;
for (i=0;i<MAX_GEN;i++)
    {
    if (char_table[i]==character) return (i);
    }
return(-1);
}


/*-----------------------------------------------------------------------------
   Prozedur, die einen Relator normalisiert, d.h. sie macht ihn lexikalisch
   so klein wie moeglich. Diese Version behandelt nur einen Relator.
   Parameter:
      relators_string - Relator.
      relators_length - Laenge des Relators.
      relatornum - Anzahl Relatoren
      char_table - Charactertabelle
----------------------------------------------------------------------------*/
void normalize_single ( relators_string, relators_length, char_table)
char* relators_string;
int relators_length;
char* char_table;
    {
    int block_start;
    int first_block;
    int first_block_start;
    int first_block_val;
    int max_block;
    int max_block_start;
    int max_block_val;
    int block_val;
    int position;
    int j;
    int l;
    int is_smaller;
    block_val=0;
    block_start=0;
    first_block=0;
    first_block_start=0;
    first_block_val = 100;
    max_block=0;
    max_block_start=0;
    max_block_val = 100;
    for (j=0;j<relators_length;j++)
        {
        if (index_table(char_table,relators_string[j])!=block_val)
           {
           if (((j-block_start)>=max_block) && (max_block_val>=block_val))
              {
              if (first_block==0) 
                 {
                 first_block = j-block_start;
                 first_block_val = block_val;
                 first_block_start = block_start;
                 }
              /* max block wird nur neu gesetzt (bei Gleichheit) ,wenn
                 der String hinter dem Block lexikalisch kleiner ist als
                 der alte */
              if (((j-block_start)>max_block) || (max_block_val>block_val))
                 {
                 max_block = j-block_start;
                 max_block_val = block_val;
                 max_block_start = block_start;
                 }
              else
              /* gefundene Bloecke sind gleich */
                 {
                 is_smaller = 0;
                 for (l=0;(l<relators_length)&&(!is_smaller);l++)
                     {
                     if (relators_string[(l+max_block_start)%relators_length] != relators_string[(l+block_start)%relators_length]) 
                        {
                        if (index_table(char_table,relators_string
                           [(l+block_start)%relators_length])
                           <index_table(char_table,relators_string
                            [(l+max_block_start)%relators_length]))
                           {
                           is_smaller = 1;
                           }
                        }
                     }
                 if (is_smaller)
                    {
                    max_block = j-block_start;
                    max_block_val = block_val;
                    max_block_start = block_start;
                    }
                 }
              }
           block_start = j;
           block_val = index_table(char_table, relators_string[j]);
           } 
        } 
        if (first_block_val==block_val)
           {
           if ((first_block+block_start-relators_length>max_block)
              && (max_block_val>=block_val))
              {
              if (((first_block+block_start-relators_length)>max_block)
                 || (max_block_val>block_val))
                 {
                 max_block = first_block+block_start-relators_length;
                 max_block_start = block_start;
                 max_block_val = block_val;
                 }
              else
              /* gefundene Bloecke sind gleich */
                 {
                 is_smaller = 0;
                 for (l=0;(l<relators_length)&&(!is_smaller);l++)
                     {
                     if (relators_string[(l+max_block_start)%relators_length] != relators_string[(l+block_start)%relators_length]) 
                        {
                        if (index_table(char_table,relators_string
                           [(l+block_start)%relators_length])
                           <index_table(char_table,relators_string
                            [(l+max_block_start)%relators_length]))
                           {
                           is_smaller = 1;
                           }
                        }
                     }
                 if (is_smaller)
                    {
                    max_block = j-block_start;
                    max_block_val = block_val;
                    max_block_start = block_start;
                    }
                 }
              }
           }
        shift_cyclic_left(relators_string,relators_length,
          max_block_start);
    }


/*-----------------------------------------------------------------------------
   Prozedur, die einen Relator normalisiert, d.h. sie macht ihn lexikalisch
   so klein wie moeglich. Diese Version behandelt alle Relatoren.
   Diese Prozedur ist notwendig, um die Relatoren vergleichbar zu machen.
   Parameter:
      relators_string - Liste der Relatoren.
      relators_length - Liste der Laengen der Relatoren.
      relatornum - Anzahl Relatoren.
      char_table - Charactertabelle
----------------------------------------------------------------------------*/
void normalize (relators_string, relators_length, relatornum, char_table)
char** relators_string;
int* relators_length;
int relatornum;
char* char_table;
{
int block_start;
int first_block;
int first_block_start;
int first_block_val;
int max_block;
int max_block_start;
int max_block_val;
int block_val;
int position;
int i;
int j;
int l;
int is_smaller;
for (i=0;i<relatornum;i++)
    {
    block_val=0;
    block_start=0;
    first_block=0;
    first_block_start=0;
    first_block_val = 100;
    max_block=0;
    max_block_start=0;
    max_block_val = 100;
    for (j=0;j<relators_length[i];j++)
        {
        if (index_table(char_table,relators_string[i][j])!=block_val)
           {
           if (((j-block_start)>=max_block) && (max_block_val>=block_val))
              {
              if (first_block==0) 
                 {
                 first_block = j-block_start;
                 first_block_val = block_val;
                 first_block_start = block_start;
                 }
              /* max block wird nur neu gesetzt (bei Gleichheit) ,wenn
                 der String hinter dem Block lexikalisch kleiner ist als
                 der alte */
              if (((j-block_start)>max_block) || (max_block_val>block_val))
                 {
                 max_block = j-block_start;
                 max_block_val = block_val;
                 max_block_start = block_start;
                 }
              else
              /* gefundene Bloecke sind gleich */
                 {
                 is_smaller = 0;
                 for (l=0;(l<relators_length[i])&&(!is_smaller);l++)
                     {
                     if (relators_string[i][(l+max_block_start)%relators_length[i]] != relators_string[i][(l+block_start)%relators_length[i]]) 
                        {
                        if (index_table(char_table,relators_string[i]
                           [(l+block_start)%relators_length[i]])
                           <index_table(char_table,relators_string[i]
                            [(l+max_block_start)%relators_length[i]]))
                           {
                           is_smaller = 1;
                           }
                        }
                     }
                 if (is_smaller)
                    {
                    max_block = j-block_start;
                    max_block_val = block_val;
                    max_block_start = block_start;
                    }
                 }
              }
           block_start = j;
           block_val = index_table(char_table, relators_string[i][j]);
           } 
        } 
        if (first_block_val==block_val)
           {
           if ((first_block+block_start-relators_length[i]>max_block)
              && (max_block_val>=block_val))
              {
              if (((first_block+block_start-relators_length[i])>max_block)
                 || (max_block_val>block_val))
                 {
                 max_block = first_block+block_start-relators_length[i];
                 max_block_start = block_start;
                 max_block_val = block_val;
                 }
              else
              /* gefundene Bloecke sind gleich */
                 {
                 is_smaller = 0;
                 for (l=0;(l<relators_length[i])&&(!is_smaller);l++)
                     {
                     if (relators_string[i][(l+max_block_start)%relators_length[i]] != relators_string[i][(l+block_start)%relators_length[i]]) 
                        {
                        if (index_table(char_table,relators_string[i]
                           [(l+block_start)%relators_length[i]])
                           <index_table(char_table,relators_string[i]
                            [(l+max_block_start)%relators_length[i]]))
                           {
                           is_smaller = 1;
                           }
                        }
                     }
                 if (is_smaller)
                    {
                    max_block = j-block_start;
                    max_block_val = block_val;
                    max_block_start = block_start;
                    }
                 }
              }
           }
        shift_cyclic_left(relators_string[i],relators_length[i], max_block_start);
    }
}
 
/*---------------------------------------------------------------------------
   Funktion zur Berechnung des Inverses eines Relators.
   Parameter:
     relators_string - Relator
     relators_length - Relatorlaenge
     orders - Ordnungen der Generatoren.
     char_table - Zuordnungstabelle Character - Indizes. 
   Rueckgabewert:
     Zeiger auf den String, der das Inverse beinhaltet.
---------------------------------------------------------------------------*/
char* inverse_relator (relators_string, relators_length,
                     orders, char_table)
char* relators_string;
int relators_length;
int* orders;
char* char_table;
{
int i;
int j;
char* invers;
char block_char;
int block_count;
int invers_len;
invers = (char*) malloc(sizeof(char)*relators_length*orders[0]);
invers_len = 0;
invers[0] = '\0';
block_char = relators_string[relators_length-1];
block_count = 1;
for (i=relators_length-2;i>=0;i--)
    {
    if (relators_string[i] == block_char)
       {
       block_count++;
       }
    else
       {
       for (j=0;j<(block_count*(orders[index_table(char_table, block_char)]-1))
              %orders[index_table(char_table, block_char)];
              j++)
           {
           invers[invers_len] = relators_string[i+1]; 
           invers_len++;
           if (invers_len>=orders[0]*relators_length)
              {
               invers = (char*) realloc (invers, sizeof(char)*(invers_len+1));
              }
           }
       block_count = 1;
       block_char = relators_string[i];
       invers[invers_len] = '\0';
       }
    }
if (block_count!=relators_length)
   {
   for (j=0;j<(block_count*(orders[index_table(char_table, block_char)]-1))
          %orders[index_table(char_table, block_char)];
          j++)
      {
      invers[invers_len] = relators_string[0]; 
      invers_len++;
      if (invers_len>=orders[0]*relators_length)
         {
          invers = (char*) realloc (invers, sizeof(char)*(invers_len+1));
         }
      }
   invers[invers_len] = '\0';
   }
else
   {
   for (j=0;j<block_count;j++)
      {
      invers[invers_len] = relators_string[0]; 
      invers_len++;
      if (invers_len>orders[0]*relators_length)
         {
          invers = (char*) realloc (invers, sizeof(char)*(invers_len+1));
         }
      }
   invers[invers_len] = '\0';
   }
return(invers);
}

/*---------------------------------------------------------------------------
   Funktion zum Loeschen der Relatoren, deren Inverses bereits in der 
   Praesentation liegt.
   Parameter:
   Rueckgabewert:
      Neue Anzahl Relatoren.
---------------------------------------------------------------------------*/
int eliminate_inverse( relators_string, relators_length, relator_count,
                       orders, char_table)
char** relators_string;
int* relators_length;
int relator_count;
int* orders;
char* char_table;
{
int i;
int j;
int k;
int is_equal;
int ende;
char* inverse_string;
int inverse_string_length;
for (i=0;i<relator_count;i++)
    {
    inverse_string = inverse_relator (relators_string[i], relators_length[i],
                     orders, char_table);
    inverse_string_length = strlen(inverse_string);
    normalize_single(inverse_string,inverse_string_length, char_table); 
    ende =0;
    for (j=i+1;(j<relator_count)&&(!ende);j++)
        {
        if (inverse_string_length==relators_length[j])
           {
           is_equal = 1;
           for (k=0;k<relators_length[j];k++)
               {
               if (relators_string[j][k]!=inverse_string[k])
                  {
                  is_equal = 0;
                  }
               }
           if (is_equal)
              {
              ende = 1;
              if (relators_length[j]<relators_length[i])
                 {
                 delete_relator_s (relators_string, relators_length, 
                       relator_count, i);
                 }
              else
                 {
                 delete_relator_s (relators_string, relators_length, 
                       relator_count, j);
                 }
              relator_count--;
              }
           }
        }
    }
return(relator_count);
}

/*----------------------------------------------------------------------------
   Prozedur, die die Laenge einer Praesentation verkuerzen soll.
   Parameter:
        presentation - In- und Output; Beinhaltet die Praesentation, die 
                       verkuerzt werden soll.
        orders - Liste der Ordnungen der Generatoren der Gruppe (in der 
                 gleichen Reihenfolge, wie in der Liste 
                 presentation->generators.
----------------------------------------------------------------------------*/
void shorten_presentation(presentation,  orders)
presentation_TYP* presentation;
int* orders;
{
int gen_no;
int j;
int* relators_length;
int i;
char** relators_string;
char* char_table;
int relator_count;
gen_no = presentation->generators->firstfree;
/* Definition der Charactertabelle */
char_table= (char*) malloc ((MAX_GEN+1)*sizeof(char));
char_table[0] = '0';
char_table[1] = '1';
char_table[2] = '2';
char_table[3] = '3';
char_table[4] = '4';
char_table[5] = '5';
char_table[6] = '6';
char_table[7] = '7';
char_table[8] = '8';
char_table[9] = '9';
char_table[10] = 'a';
char_table[11] = 'b';
char_table[12] = 'c';
char_table[13] = 'd';
char_table[14] = 'e';
char_table[15] = 'f';
char_table[16] = 'g';
char_table[17] = 'h';
char_table[18] = 'i';
char_table[19] = 'j';
char_table[20] = 'k';
char_table[21] = 'l';
char_table[22] = 'm';
char_table[23] = 'n';
char_table[24] = 'p';
char_table[25] = '\0';
/* Relatoren in Strings konvertieren */
relators_length = (int*) malloc (sizeof(int)*presentation->norelators);
relators_string = (char**) malloc (sizeof(char*)*presentation->norelators);
for (i=0;i<presentation->norelators;i++)
    {
    relators_string[i] = (char*) malloc (sizeof(char)*(presentation->relators[i].lhsnproduct+1));
    relators_length[i] = presentation->relators[i].lhsnproduct;
    for (j=0;j<relators_length[i];j++)
        {
        relators_string[i][j] = char_table[presentation->relators[i].lhsproduct[j]];
        }
    relators_string[i][relators_length[i]] = '\0';
    }
/**********
for (i=0;i<presentation->norelators;i++)
    {
    printf(" relator %d: %s \n", i, relators_string[i]);
    }
***********/
/* zunaechst werden doppelte Relatoren eliminiert */
relator_count = presentation->norelators;
fprintf(stderr," no relators before eliminating duplicates: %d \n",
       relator_count);
relator_count = eliminate_duplicates( relators_string, relators_length, relator_count);
fprintf(stderr," no relators after eliminating duplicates (first pass): %d \n",
       relator_count);
/*****
for (i=0;i<relator_count;i++)
    {
    printf(" relator %d: %s \n", i, relators_string[i]);
    }
******/
/* Normalisieren der Relatoren und erneuetes Eliminieren von Duplikaten */
normalize( relators_string, relators_length, relator_count, char_table);
relator_count = eliminate_duplicates( relators_string, relators_length, relator_count);
fprintf(stderr," no relators after eliminating duplicates (second pass): %d \n",
       relator_count);
/* Loeschen der Relatoren, deren Inverses bereits in der Liste ist,
   bzw. loeschen des laengeren von beiden */
relator_count = eliminate_inverse( relators_string, relators_length, 
                 relator_count, orders, char_table);
fprintf(stderr," no relators after eliminating duplicates (inverses)");
fprintf(stderr," (third pass): %d \n",relator_count);
/* Rueckumwandeln der Characterstrings in die alte Form der Praesentation */
for (i=relator_count;i<presentation->norelators;i++)
    {
    free(presentation->relators[i].lhsproduct);
    presentation->relators[i].lhsnproduct = 0;
    }
presentation->norelators = relator_count;
for (i=0;i<relator_count;i++)
    {
    presentation->relators[i].lhsproduct = (int*) realloc (presentation->relators[i].lhsproduct, sizeof(int)*relators_length[i]);
    presentation->relators[i].lhsnproduct = relators_length[i];
    for (j=0;j<relators_length[i];j++)
        {
        presentation->relators[i].lhsproduct[j] = index_table(char_table,relators_string[i][j]); 
        }
    }
/*
for (i=0;i<presentation->norelators;i++)
    {
    free (relators_string[i]);
    }
*/
free (relators_string);
free (relators_length);
}

