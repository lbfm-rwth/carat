/* some definitions of constants */

/* the maximal entry in the short vectors should not exceed MAXENTRY to
   avoid overflow */
#define	MAXENTRY	10000
/* it is checked that the product of the maximal entry in the short vectors
   and the maximal entry in the products of Gram-matrices with short vectors
   does not exceed MAXNORM, since the program works just with int's */
#define	MAXNORM		100000000
#define	STRLEN		80
/* since the LLL-reduction works with a floating-point model, there might be
   some rounding error, which should not reach EPS */
#define	EPS		0.001
/* the constant for the LLL-condition, which should be < 1, but 0.75 is the
   standard value, a higher value might yield a slightly better result but
   increases the running time */
#define	LLL_CONST	0.75
#define DEF_PRIME	16001


/* structure to hold the generators of the group according to the stabilizer
   chain */
typedef struct {
	int	dim;
	int	*ord;
	int	*ng;
	int	*nsg;
	int	****g;
	} group;

/* structure for a list of vectors */
typedef struct {
	int	dim;
	int	len;
	int	prime;
	int	n;
	int	**v;
	} veclist;

/* structure for the invariant forms and their products with the list of
   short vectors */
typedef struct {
	int	dim;
	int	n;
	int	***A;
	int	***v;
	} invar;

/* structure for the diagonal of the fingerprint, the new order and the indices
   of the standard-base in the list of short vectors */
typedef struct {
	int	*diag;
	int	*per;
	int	*e;
	} fpstruct;

/* structure for the scalar product combinations, the transformation matrices
   between the vector sums and that lattice bases and the scalar products of
   the lattice bases */
typedef struct {
	veclist	list;
	int	rank;
	int	**trans;
	int	**coef;
	int	***F;
	} scpcomb;

/* structure for a Bacher-polynomial */
typedef struct {
	int	mind;
	int	maxd;
	int	sum;
	int	*coef;
	} bachpol;

/* structure for the option flags */
typedef struct {
	int	DEPTH;
	int	STAB;
	int	BACH[3];
	int	BACHDEP;
	int	BACHSCP;
	int	GEN;
	int	PRINT;
	} flagstruct;

/* functions in auttools.c */
static	int	cand();
static	void	autom();
static	int	aut();

/* functions in bachtools.c */
static	void	bacher();
static	int	bachcomp();
static	void	fputbach();

/* functions in iotools.c */
static  void	getflags();
static	bravais_TYP 	*putgens();
static	void	putord();
static	matrix_TYP	*putiso();

/* functions in isotools.c */
static	int	isocand();
static	matrix_TYP 	*bs_isometry();
static	int	iso();
static	int	isostab();

/* functions in lattools.c */
static	int	lll();
static	void	initialize();
static	int	red();
static	void	check();
static	void	decrease();
static	void	interchange();
static	int	iround();

/* functions in mattools.c */
static	void	vecmatmul();
static	void	matmul();
static	int	scp();
static	int	sscp();
static	void	psolve();
static	void	pgauss();
static	int	isprime();

/* functions in orbtools.c */
static int	operate();
static int	orbit();
static int	orbitlen();
static int	delete();
static void	stab();
static void	matgen();
static void	stabil();

/* functions in preproc.c */
static void	checkvecs();
static int	checkgen();
static void	fingerprint();
static int	possible();
static void	scpvector();
static void	scpvecs();
static void	base();
static void	coef();
static void	scpforms();

/* functions in sorttools.c */
static	int	comp();
static	int	numberof();
static	void	sortvecs();
static	void	quicksort();

/* functions in perfecttools.c */
static int normal_aut_test();
