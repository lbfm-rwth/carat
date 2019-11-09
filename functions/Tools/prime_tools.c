#include "typedef.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: prime_tools.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*
@-------------------------------------------------------------------------
@ This file exports the global variables
@ int act_prime  - This the actual used prime_number
@
@ int (*S)(int,int), (*P)(int,int)
@ S is the function for addition modulo act_prime and
@ P the function for multiplication
@
@ void cleanup_prime();
@             releases all memory that was allocated by init_prime(). Only 
@             useful for debugging memory allocation.
@ void init_prime( int prime);
@    This function MUST be called before the functions S(a,b) and P(a,b)
@    can be used, since otherwise they are not defined.
@    The argument of init_prime is the new prime number, modulo which
@    the calculations shoud be done, i.e. act_prime is set to prime.
@
@-------------------------------------------------------------------------
| Dieses file export die globalen Variablen
|
| int act_prime;
|               Dies ist die gerade benutzte Primzahl
|
| int (*S)(int, int), (*P)(int,int);
|               Dies sind die Multiplikationsfunktionen fuer das Rechenen
|               modulo act_prime;
|
| Es werden zwei Funktionen definiert, naemlich
|
| void init_prime( int prime );
|              Diese Funktion \underline{MUSS} aufgerufen werden, bevor die
|              Rechenfunktionen S(a,b) und P(a,b) benutzt werden koennen,
|              da diese sonst nicht definiert sind. Argument von
|              init_prime() ist die neue Primzahl, modulo derer gerechnet
|              werden soll, das heisst, act_prime wird auf prime gesetzt
|
| void cleanup_prime();
|             releases all memory that was allocated by init_prime(). Only 
|             useful for debugging memory allocation.
|
*/
#include "tools.h"

/*
 *  local variables
 */
#define SHORTPRIME_LIMIT 300 /* Grenze fuer Rechnen mit Tabelle  */

typedef struct {
  int s_number, *s_primes; /* Anzahl kleiner Pz; Liste der ... */
  int l_number, *l_primes; /* Anzahl grosser Pz; Liste der ... */
  int ***ADD, ***MUL;       /* Rechen-Tabellen fuer kleine Pz.  */
  int **LOG, **EXP;         /* p-Logarithmen fuer grosse Pz.    */
  } table;

static int dummy( int a, int b);

static int **P_ADD = NULL, **P_MUL = NULL;               /* aktuelle Rechen_Tabelle */
static int *P_LOG = NULL, *P_EXP = NULL;                 /* aktuelle p-Logarithmen  */
static table *prime_table = NULL;

/* 
 * global variables
 */
int (*S)(int, int)= dummy, (*P)(int, int)= dummy;        /* aktuelle Funktion       */
int act_prime= -1;

/*}}}  */
/*{{{  dummy*/
static int 
dummy (int a, int b) 
{ 
  fprintf(stderr,"You tried to use GF(p)-Arithmetic without calling init_prime(act_prime).\n");
  fprintf(stderr,"Currently, act_prime is set to %d, which is apparently not prime.\n", act_prime);
  exit(3);
  return 0;
}

/*}}}  */
/*{{{  short_P*/
static int 
short_P (int i, int j)
{
  return(P_MUL[i][j]);
}

/*}}}  */
/*{{{  short_S*/
static int 
short_S (int i, int j)
{
  return(P_ADD[i][j]);
}

/*}}}  */
/*{{{  int_P*/
static int 
int_P (int i, int j)
{
  return((i*j ==0 ) ? 0 : P_EXP[(P_LOG[i] + P_LOG[j]) % (act_prime-1)]);
}

/*}}}  */
/*{{{  int_S*/
static int 
int_S (int i, int j)
{
/* die "5" ist nur ein Sicherheitsfaktor */
  return((i+j+5*act_prime)%act_prime);
}

/*}}}  */
/*{{{  init_short_prime*/
static void 
init_short_prime (void)
{
int j,k;

  P_ADD = (int **)malloc((unsigned)act_prime*sizeof(int *));
  P_MUL = (int **)malloc((unsigned)act_prime*sizeof(int *));
  for ( j = 0; j < act_prime; j ++) {
    P_MUL[j] = (int *)malloc((unsigned)(2*act_prime-1)*sizeof(int));
    P_ADD[j] = (int *)malloc((unsigned)(2*act_prime-1)*sizeof(int));
    P_MUL[j] += act_prime-1;
    P_ADD[j] += act_prime-1;
  }
  for ( j = 0; j < act_prime; j ++) {
    for ( k = 0; k <= j; k++) {
      P_MUL[j][k] = P_MUL[k][j] = (j * k) % act_prime;
      if ( j != 0 ) {
        P_MUL[P_MUL[j][k]][-j] = k;
      }
      if ( k != 0 ) {
        P_MUL[P_MUL[j][k]][-k] = j;
      }
      if((k != 0) && (P_MUL[j][k] == 0)) {
        fprintf(stderr,"init_short_prime: Error Non-prime (%d) in makefield occured\n", act_prime);
      }
      P_ADD[j][k] = P_ADD[k][j] = (j + k) % act_prime;
      P_ADD[P_ADD[j][k]][-j] = k;
      P_ADD[P_ADD[j][k]][-k] = j;
    }
  }
}
/*}}}  */
/*{{{  init_int_prime*/
static void 
init_int_prime (void)
{
int i,j,k,gt;
int  *p_tmp;
int flag;

  P_LOG = (int *)calloc(2*act_prime+1,sizeof(int));
  P_EXP = (int *)calloc(act_prime,sizeof(int));
  p_tmp = (int *)calloc(act_prime+1,sizeof(int));
  P_LOG += act_prime;
  
  /* find a generator of the cyclic multiplicative group of Fp */
  
  for(i = 1; i <= act_prime; i++) {
    p_tmp[i] = 1;
  }
  
  k = gt = (act_prime-1)/2;
  p_tmp[k] = 0;
  j = 1;
  flag = TRUE;
  while(flag) {
    do {
      if((k = (k * gt)%act_prime) == 0) {
        fprintf(stderr,"Error: Nullteiler\n");
        exit(3);
      }
      p_tmp[k] = 0;
      j++;
    } while(k != 1);
    flag = (j == act_prime-1) ? FALSE : TRUE;
    if(flag) {
      j = 1;
      for(i = 1; i <=act_prime; i++) {
        if (p_tmp[i] == 1) {
          gt = i;
          break;
        }
      }
      if(i > act_prime) {
      fprintf(stderr,"No generator for multiplicative group found");
        exit(3);
      }
      k = gt;
      p_tmp[k] = 0;
    }
  }
  
  k = j = 1;
  do {
    k = (k * gt)%act_prime;
    P_LOG[k] = j;
    P_EXP[j] = k;
    P_LOG[-k] = act_prime-j-1;
    j++;
  } while(k != 1);
  P_LOG[1] = 0;
  P_EXP[0] = 1;
  free(p_tmp);
}
/*}}}  */
/*{{{  cleanup_prime*/
void 
cleanup_prime (void)
{                    
int i, j;

  if ( prime_table != NULL ) {
    for ( i=1; i < prime_table->s_number; i ++ ) {
      act_prime = prime_table->s_primes[i];
      for ( j = 0; j < act_prime; j ++) {
        prime_table->ADD[i][j] -= act_prime - 1;
        prime_table->MUL[i][j] -= act_prime - 1;
        free( prime_table->ADD[i][j] );
        free( prime_table->MUL[i][j] );
      }                         
      free( prime_table->ADD[i] );
      free( prime_table->MUL[i] );
    }
    free( prime_table->s_primes );
    free( prime_table->ADD );
    free( prime_table->MUL );
    
    for ( i=1; i < prime_table->l_number; i ++ ) {
      act_prime = prime_table->l_primes[i];
      prime_table->LOG[i] -= act_prime;
      free( prime_table->LOG[i] );
      free( prime_table->EXP[i] );
    }
    free( prime_table->l_primes );
    free( prime_table->LOG );
    free( prime_table->EXP );
    free( prime_table );
  }
  act_prime = -1;
  prime_table = NULL;
  P_LOG = NULL;
  P_EXP = NULL;
  P_ADD = NULL;
  P_MUL = NULL;


}
/*}}}  */
/*{{{  init_prime*/
void 
init_prime (int prime)
{
int j, number;

  if (act_prime == -1) {
    prime_table = (table *)malloc(sizeof(table));
    prime_table->s_number = 1;
    prime_table->l_number = 1;
    prime_table->s_primes= (int   *)calloc(1,sizeof(int));
    prime_table->l_primes= (int   *)calloc(1,sizeof(int));
    prime_table->ADD   = (int ***)malloc(1*sizeof(int **));
    prime_table->MUL   = (int ***)malloc(1*sizeof(int **));
    prime_table->LOG   = (int  **)malloc(1*sizeof(int  *));
    prime_table->EXP   = (int  **)malloc(1*sizeof(int  *));
  }
  j = 0;
  if (prime < SHORTPRIME_LIMIT) {
    number = prime_table->s_number;
    while (j < number) {
      if (prime_table->s_primes[j] == prime) {
        P_ADD = prime_table->ADD[j];
        P_MUL = prime_table->MUL [j];
        act_prime = prime;
        return;
      }
      j++;
    }
    number++;
    prime_table->s_primes= (int *)realloc(prime_table->s_primes,
                      number*sizeof(int));
    prime_table->ADD = (int ***)realloc(prime_table->ADD,
                      number*sizeof(int **));
    prime_table->MUL = (int ***)realloc(prime_table->MUL ,
                      number*sizeof(int **));
    prime_table->s_primes[number-1] = act_prime = prime;
    init_short_prime();
    prime_table->ADD[number-1] = P_ADD;
    prime_table->MUL[number-1] = P_MUL;
    prime_table->s_number = number;
    S = (int (*)(int, int))(short_S);
    P = (int (*)(int, int))(short_P);
  } else {
    number = prime_table->l_number;
    while (j < number) {
      if (prime_table->l_primes[j] == prime) {
        P_LOG = prime_table->LOG[j];
        P_EXP = prime_table->EXP[j];
        act_prime = prime;
        return;
      }
      j++;
    }
    number++;
    prime_table->l_primes= (int   *)realloc(prime_table->l_primes,
                      number*sizeof(int));
    prime_table->LOG = (int **)realloc(prime_table->LOG,
                      number*sizeof(int *));
    prime_table->EXP = (int **)realloc(prime_table->EXP,
                      number*sizeof(int *));
    prime_table->l_primes[number-1] = act_prime = prime;
    init_int_prime();
    prime_table->LOG[number-1] = P_LOG;
    prime_table->EXP[number-1] = P_EXP;
    prime_table->l_number = number;
    S = (int (*)(int, int))(int_S);
    P = (int (*)(int, int))(int_P);
  }     


  return;
}

/*}}}  */
