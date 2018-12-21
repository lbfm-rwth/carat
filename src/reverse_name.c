#include "ZZ.h"
#include "typedef.h"
#include "getput.h"
#include "name.h"
#include "bravais.h"
#include "datei.h"
#include "matrix.h"
#include "voronoi.h"
#include "autgrp.h"
#include "symm.h"
#include "base.h"
#include "zass.h"
#include "gmp.h"
#include "longtools.h"

int SFLAG;
int INFO_LEVEL;
boolean GRAPH = FALSE;


int main (int argc, char *argv[])
{

  bravais_TYP *R;

  char qname[1024];


  int zname[2],
    i, c;

  MP_INT aff_name;

  char comment[1024];

  char *affstring;

  
  
  read_header (argc, argv);


  if (is_option('h'))
    INFO_LEVEL = optionnumber('h');

  if (INFO_LEVEL == 8)
    SFLAG = 1;

  if ((is_option('h') && INFO_LEVEL != 8) || FILEANZ > 0)
    {
      printf("Usage: %s [-h] [-c] [-i]\n",argv[0]);
      printf("\n");
      printf("The program enables you to construct a space group corresponding\n");
      printf("to a name as given by CARAT. The program will read the various\n");
      printf("components of the name for the desired space group from stdin.\n");
      printf("It will then output generators for the space group, which are\n");
      printf("understood to generate the space group together with Z^n.\n");
      printf("Please note: not every name is a valid name!\n");
      printf("\n");
      printf("Options:\n");
      printf("-h     : gives you this help.\n");
      printf("-c     : do not check the name given, ie. verify that it is a valid name.\n");
      printf("         WARNING: this could lead to a wrong name in the header of the\n");
      printf("         resulting group\n");
      printf("-i     : ignore that the name is invalid, and give a space group at least\n");
      printf("         in the desired Q-class, and if exists one in the desired Z-class.\n");
      printf("         If given without -c the resulting group will also indicate the valid\n");
      printf("         name.\n");
      printf("Note: The Q-classes in               correspond to the Q-classes in degree\n");
      printf("   min.1 max.1                                                     1\n");
      printf("   min.2-5,group.1-4, max.2-3                                      2\n");
      printf("   min.6-14,group.5-25,max.4-5                                     3\n");
      printf("   min.15-57,group.26-205,max.6-9                                  4\n");
      printf("   min.58-169,group.206-1042,max.10-15                             5\n");
      printf("   min.170-667,group.1043-7636,max.16-27                           6\n");
      printf("\n");
      printf("Cf. Name, Q_catalog, QtoZ, Extensions.\n");

      if (FILEANZ == 0)
         exit(0);
      else
         exit(31);
    }


  fprintf(stderr,"qname:\n");
  c=scanf("%s",qname);

  fprintf(stderr,"zname: \n");
  c=scanf("%d %d",zname,zname+1);

  fprintf(stderr,"affname: \n");
  mpz_init(&aff_name);
  mpz_inp_str(&aff_name,stdin,10);


  if (is_option('c')){
     if (is_option('i')){
        i = 2;
     }
     else{
        i = 0;
     }
  }
  else{
     if (is_option('i')){
        i = 3;
     }
     else{
        i = 1;
     }
  }

  R = reverse_name(qname, zname, aff_name, i, is_option('i'), &affstring);
  sprintf(comment,"standard group with name %s %d %d %s",
                   qname,zname[0],zname[1],affstring);
  put_bravais(R,NULL,comment);




  free(affstring);
  free_bravais(R);
  if (INFO_LEVEL == 8) pointer_statistics(0,0);

  exit(0);
}
