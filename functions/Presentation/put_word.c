#include "typedef.h"

/***************************************************************************
@
@---------------------------------------------------------------------------
@ FILE: put_word.c
@---------------------------------------------------------------------------
@
****************************************************************************/

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ void put_word(int *w,
@               char *O)
@ outputs the word w to stdout.
@ the string O contatins the options. If O contains a 'G', the output
@ is done in gap format, if O contains a 'M', the output is
@ done in matrix format. Both is possible simultanously.
@---------------------------------------------------------------------------
@
****************************************************************************/
void put_word(int *w,
              char *O)
{

  int i;

  if (strchr(O,'G')){
    if (w[0] == 0){
      printf("g[1]*g[1]^(-1)");
    }
    else{
      for(i=1;i<w[0];i++){
        if (w[i] > 0){
	  printf("g[%d]*",w[i]);
	}
        else{
	  printf("g[%d]^(-1)*",-w[i]);
	}
        if (!(i % 10)) printf("\n");
      }
      if (w[i] > 0){
	printf("g[%d]",w[i]);
      }
      else{
	printf("g[%d]^(-1)",-w[i]);
      }
    }

    printf("\n");
  }

  if (strchr(O,'M')){
    if (w[0] == 0){
      printf("c 1 -1");
    }
    else{
      printf("c ");
      for(i=1;i<=w[0];i++)
	printf("%d ",w[i]);
    }

    printf("\n");
  }


  fflush(stdout);

  return;

} /* put_word(...) */

