#ifdef __cplusplus
extern "C" {
#endif


#ifndef _TOOLS_H_
#define _TOOLS_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: carat_exit.c 
\*-------------------------------------------------------------*/
extern void carat_exit(char *fmt, ...);

/*-------------------------------------------------------------*\
| FILE: chin_remainder.c 
\*-------------------------------------------------------------*/
extern int chin_remainder(int x1, int x2, int p1, int p2);

/*-------------------------------------------------------------*\
| FILE: intpow.c 
\*-------------------------------------------------------------*/
extern int intpow(int a,int b);

/*-------------------------------------------------------------*\
| FILE: itoasc.c 
\*-------------------------------------------------------------*/
extern void itoasc(int n, char s[]);

/*-------------------------------------------------------------*\
| FILE: malloc2dim.c 
\*-------------------------------------------------------------*/
extern char **calloc2dim(int r,int c,int size);
extern char **malloc2dim(int r,int c,int size);
extern void memcpy2dim(char **dest, char **src, int r,int c,int size);
extern void memset2dim(char **dest, int r, int c, int size, char *value);
extern void free2dim(char **old, int rows);

/*-------------------------------------------------------------*\
| FILE: mindiv.c 
\*-------------------------------------------------------------*/
extern int min_div( int a, int b) ;

/*-------------------------------------------------------------*\
| FILE: ovfl_mul.c 
\*-------------------------------------------------------------*/
extern int ovfl_mul( int a, int b);

/*-------------------------------------------------------------*\
| FILE: prime_tools.c 
\*-------------------------------------------------------------*/
extern int act_prime;

extern int (*S)(int, int);
extern int (*P)(int, int);
extern void cleanup_prime();
extern void init_prime ( int prime);

/*-------------------------------------------------------------*\
| FILE: ramdom.c
\*-------------------------------------------------------------*/
extern int random_own();

/*-------------------------------------------------------------*\
| FILE: tools.c 
\*-------------------------------------------------------------*/
extern rational Zero;
extern rational One;

extern int GGT (int _a, int _b);
extern int KGV(  int a,  int b );
extern void Normal (rational *a);
extern void Normal2 ( int *_z, int *_n );
extern void rat_add( int *az, int *an, int bz, int bn );
extern int *factorize_new( int zahl, int *erg);
extern int *factorize( int zahl);
extern void gcd_darstell(int a1, int a2, int *v1, int *v2, int *gcd);
extern int p_inv(int a, int p);
extern int signum(int a);

#else
/*-------------------------------------------------------------*\
| FILE: carat_exit.c 
\*-------------------------------------------------------------*/
extern void carat_exit();

/*-------------------------------------------------------------*\
| FILE: chin_remainder.c 
\*-------------------------------------------------------------*/
extern int chin_remainder();

/*-------------------------------------------------------------*\
| FILE: intpow.c 
\*-------------------------------------------------------------*/
extern int intpow();

/*-------------------------------------------------------------*\
| FILE: itoasc.c 
\*-------------------------------------------------------------*/
extern void itoasc();

/*-------------------------------------------------------------*\
| FILE: malloc2dim.c 
\*-------------------------------------------------------------*/
extern char **calloc2dim();
extern char **malloc2dim();
extern void memcpy2dim();
extern void memset2dim();
extern void free2dim();

/*-------------------------------------------------------------*\
| FILE: mindiv.c 
\*-------------------------------------------------------------*/
extern int min_div();

/*-------------------------------------------------------------*\
| FILE: ovfl_mul.c 
\*-------------------------------------------------------------*/
extern int ovfl_mul();

/*-------------------------------------------------------------*\
| FILE: prime_tools.c 
\*-------------------------------------------------------------*/
extern int act_prime;

extern int (*S)();
extern int (*P)();
extern void cleanup_prime();
extern void init_prime ();

/*-------------------------------------------------------------*\
| FILE: ramdom.c
\*-------------------------------------------------------------*/
extern int random_own();

/*-------------------------------------------------------------*\
| FILE: tools.c 
\*-------------------------------------------------------------*/
extern rational Zero;
extern rational One;

extern int GGT ();
extern int KGV();
extern void Normal ();
extern void Normal2 ();
extern void rat_add();
extern int *factorize_new();
extern int *factorize();
extern void gcd_darstell();
extern int p_inv();
extern int signum();

#endif
#endif


#ifdef __cplusplus
}
#endif


