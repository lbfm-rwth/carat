#include "typedef.h"
#include "tools.h"

int INFO_LEVEL = 0;

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: tools.c
@ containes some general small tools.
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
@
@ Exports the globale variables Zero and One.
@        rational Zero=  { 0,1 };
@        rational One =  { 1,1 };
@
@ int GGT( int a, int b );
@       calculates the greatest common divisor of a and b.
@       The result is always > 0
@
@ int KGV( int a, int b );
@       calculates the least commom multiple of a and b.
@       The result is always > 0
@
@ void Normal ( rational *a );
@       divides the numerator  and the denominator of a rational number
@       by their greatest common divisor
@
@ void Normal2 ( int *z, *n );
@       divides z and n by their greatest common divisor.
@       Caution: z and n are changed in Normal2
@
@ void rat_add ( int *az, int *an, int bz, int bn);
@                   cz       az   bz
@      calculates   --  :=   -- + --
@                   cn       an   bn
@       The result is stored in az and an.
@
@ int *factorize_new( int zahl,int *erg);
@
@ int *factorize ( int zahl );
@       The integer 'zahl' is factorized. If for example
@       zahl = p1^a1 *p2^a2,
@       erg[p1] is set to a1 and erg[p2] is set to a2,
@       all other entries of erg are 0.
@       Caution: the result is a pointer to an integer of length 100.
@                So factorize works only for intergers that
@                have prime factors smaller then 100.
@
@  void gcd_darstell( int a1, int a2, int *v1, int *v2, int *gcd);
@      calculates a presentation of the greatest commom divisor of a1 and a2:
@      gcd = v1*a1 + v2*a2
@
@  int p_inv(int a, int p)
@      calculates number ai such that a * ai is kongruent 1 modulo p (if exists)
@
@-------------------------------------------------------------------------
\**************************************************************************/

/**********************************************************************\
|
| tools.c -- allgemeine kleinere tools:
|
| int GGT( int a, int b );
| int KGV( int a, int b );
| void Normal ( rational *a );
| void Normal2 ( int *z, *n );
| void rat_add ( int *az, int *an, int bz, int bn);
| int *factorize_new( int zahl,int *erg);
| int *factorize ( int zahl );
| void gcd_darstell( int a1, int a2, int *v1, int *v2, int *gcd);
| int p_inv(int a, int p);
| int signum(int a);
|
| exportiert die globale Variable Zero. Eine Moeglichkeit, alles zum
| Absturz zu bringen, ist also z.B. die Anweisung Zero.n = 0 ...  :(
|
\**********************************************************************/

rational Zero=  { 0,1 };
rational One =  { 1,1 };

/*{{{}}}*/
/*{{{  GGT*/
int GGT (int _a, int _b)
{        
register int a= _a;
register int b= _b;
register int c;

  if ( b == 0 ) {
    return abs(a);
  } else if ( b == 1 || a == 1 ) {
    return 1;
  } else {
    while ( (c = a%b) != 0 ) {
      a = b;
      b = c;
    }
  }
  return abs(b);
}

/*}}}  */
/*{{{  KGV*/
/*
|
| int KGV( int a, int b);
|
| berechnet das kgv von a und b. Das Ergebnis ist immer >= 0!
|
 */
int 
KGV (int a, int b)
{ 
int kgv, ggt;

  a = abs(a);
  b = abs(b);
  ggt= GGT( a, b );
  if ( ggt != 0 ) 
  {
    if ( a > b )
      kgv = (a / ggt) * b; /* die Reihenfolge ist wesentlich!    */
                           /* (a*b)/ggt waere grosser Bloedsinn! */
    else
      kgv = (b / ggt) * a; /* die Reihenfolge ist wesentlich! */
  }                                                              
  else
    kgv = 0;
  return kgv;
}

/*}}}  */
/*{{{  Normal*/


void 
Normal (rational *a)
{  
register int g; 
register int n= a->n;
register int z= a->z;

  if ( n == 0 ) {
    printf ("Normal: Error: divide by zero\n");
    exit (3);
  }
  if ( z == 0 ) {
    a->n= 1;       
  } else {
    g = GGT ( z, n);
    if ( g != 1 ) {
      z /= g;
      n /= g;
    }
    if ( n < 0 ) {
      a->z = -z;
      a->n = -n;
    } else {
      a->z = z;
      a->n = n;
    }
  }
}

/*}}}  */
/*{{{  Normal2*/
void 
Normal2 (int *_z, int *_n)
{          
register int z = *_z;
register int n = *_n;
register int g;

  if ( n == 0 ) {
    printf ("Normal2: Error: divide by zero\n");
    exit (3);
  }
  if ( z == 0 ) {
    *_n = 1;
  } else {
    g = GGT ( z, n);
    if ( g != 1 ) {
      z /= g;
      n /= g;
    }
    if ( n < 0) {
      *_z = -z;
      *_n = -n;
    } else {
      *_z = z;
      *_n = n;
    }
  }
}

/*}}}  */
/*{{{  rat_add*/
/*
|
| rat_add( &a.z, &a.n, b.z, b.n );
|
| rational a, rational b: addition of a and b, result is stored in a.
|
 */
void 
rat_add (int *az, int *an, int bz, int bn)
{        
register int temp_ggt;

  /*
   *  Normal tests also for divide by zero
   */
  Normal2( az, an );
  Normal2( &bz, &bn );
  temp_ggt = GGT( *an, bn );             
  *az = (*az) * bn/temp_ggt + bz * (*an)/temp_ggt;
  if ( *an > bn ) {
    *an = (*an/temp_ggt) * bn;
  } else { 
    *an *= (bn/temp_ggt);
  }         
  Normal2( az, an );
}

/*}}}  */
/*{{{  factorize_new*/
int *
factorize_new (int zahl, int *erg)
{
int i;

  if ( erg != NULL )
  {
    for(i=0; i<100; i++)
      erg[i] = 0;
    for(i=2; i<100; i++)
    {
      while((zahl%i) == 0)
      {
        zahl = zahl/i;
        erg[i]++;
      }
    }
    if(zahl != 1)
     erg[0] = 1;
  }
  return( erg );
}

/*}}}  */
/*{{{  factorize*/
int *
factorize (int zahl)
{
int *erg;

  erg = (int *) malloc(100 *sizeof(int));
  factorize_new( zahl, erg );
  return(erg);
}

/*}}}  */
/*{{{  gcd_darstell*/
/*-----------------------------------------------------------------------*\
| gibt darstellung des ggt: gcd = v1*a1 + v2*a2                           |
\*-----------------------------------------------------------------------*/
void 
gcd_darstell (int a1, int a2, int *v1, int *v2, int *gcd)
{
 int bn,bn1,bn2,rn,rn1,rn2,q, an, an1, an2;
  rn2=a2; rn1=a1; bn2=0; bn1=1;
  an2 = 1; an1 = 0;

  if(a1 == 0 && a2 == 0)
  {
    *gcd = 0; *v1 = 0; *v2 = 0;
  }
  if(a1 == 0 && a2 != 0)
  {
    *gcd = a2; *v1 = 0; *v2 = 1;
  }
  if(a1 != 0 && a2 == 0)
  {
    *gcd = a1; *v1 = 1; *v2 = 0;
  }
  if(a1 != 0 && a2 != 0)
  {
    do
    {
       rn=rn2%rn1;
       q=(rn2-rn)/rn1;
       bn=bn2-q*bn1;
       an=an2-q*an1;
       if(rn != 0)
       {
         bn2=bn1; bn1=bn; rn2=rn1; rn1=rn; an2=an1; an1 = an;
       }
    }
    while(rn!=0);
    if(rn1 < 0)
    {
       *gcd = -rn1; *v1 = -bn1; *v2 = -an1;
    }
    if(rn1 > 0)
    {
       *gcd = rn1; *v1 = bn1; *v2 = an1;
    }
  }
}
/*}}}  */

/*-----------------------------------------------------------------------*\
| calculates number i such that a * i is kongruent 1 modulo p (if exists)
\*-----------------------------------------------------------------------*/
int 
p_inv (int a, int p)
{
   int an, an1, an2, rn, rn1, rn2, q;
   rn1 = a %p;
   if(rn1 == 0)
   {
    printf("%d is not invertible modulo %d\n", a, p);
    exit(3);
   }
   if(rn1 == 1 || rn1 == -1)
    return(rn1);
   rn2 = p; an2 = 0; an1 = 1;
   do
   {
      rn=rn2%rn1;
      q=(rn2-rn)/rn1;
      an=an2-q*an1;
      if(rn != 0)
      {
        rn2=rn1; rn1=rn; an2=an1; an1 = an;
      }
   }
   while(rn!=0);
   if(rn1  == -1)
      return( -an1);
   if(rn1  == 1)
      return( an1);
   else
   {
     printf("%d is not invertible modulo %d\n", a, p);
     exit(3);
   }
}


int 
signum (int a)
{
  if(a>0)
  return(1);
  if(a<0)
  return(-1);
 return(0);
}
