#include "typedef.h"
#include "getput.h"
#include "matrix.h"


/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: put_order.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void fput_order(outfile, divisors, ord)
@ FILE *outfile;
@ int *divisors;
@ int ord;
@ A tools to print the order of a bravais_TYP
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
fput_order (FILE *outfile, int *divisors, int ord)
{
  int tester, i, j;
  if(divisors == NULL)
  {
     if(ord <= 0)
       fprintf(outfile, "%% order of the group unknown\n");
     else
       fprintf(outfile, " = %d %% order of the group\n", ord);
  }
  else
  {
    tester = FALSE;
    if(divisors[0] != 0)
       fprintf(outfile, "%% order of the group unknown\n");
    else
    {
      for(i=0; i<100 && tester == FALSE; i++)
        if(divisors[i] != 0)
           tester = TRUE;
      if(tester == FALSE && ord == 0)
       fprintf(outfile, "%% order of the group unknown\n");
      if(tester == TRUE)
      {
        j = 0;
        for(i=2; i<100; i++)
        {
          if(divisors[i] != 0)
          {
              if(j != 0)
                fprintf(outfile, "* ");
              fprintf(outfile, "%d^%d  ", i, divisors[i]);
              j = 1;
          }
        }
          if(ord != 0)
             fprintf(outfile, " = %d ", ord);
          fprintf(outfile, "%% order of the group\n");
      }
    }
  }
}


/*{{{}}}*/
