
/* Struktur in der eine Wand mit allen notwendigen Daten abgespeichert werden
   kann. hplane beinhaltet die Gleichung der Wand. product beinhaltet die
   Darstellung des Gruppenelements als Product der Erzeuger, das auf den
   urspruenglichen Vektor angewendet wurde. nproduct gibt die groesse des
   Feldes product an. left und right sind fuer die Abspeicherung der Waende
   in einem Baum gedacht.
*/
typedef struct{
               int* hplane;
               int dist;
               int* product;
               int nproduct;
            /*    int next;*/
               int left, right;
              } wall_Typ;

/*
   Struktur fuer eine Waendeliste (als Baumstruktur (Indexverpointerung)).
   dimension gibt an welche Dimension der Raum hat in dem die Waende liegen.
   sizemult gibt an, wieviel Speicherplatz fuer die Liste bereits angefordert
   wurde (als Einheit dient hier die Konstante EXT_SIZE). firstfree ist
   der Index des ersten freien Speicherplatzes in der Liste. 
*/
typedef struct{
               wall_Typ** wall_list;
               int dimension;
               int sizemult;
               int firstfree;
            /*   int first;*/
              } wall_list_TYP;
#define nichts (-1)
