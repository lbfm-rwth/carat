#include "tools.h"

#define RANDOM_PRIME 19853
#define step 12345
int 
random_own (void)
{

  static int erg;

  erg += step;
  erg = erg % RANDOM_PRIME;

  return erg;

}
