#ifndef ARITH_TYPEDEF
#define ARITH_TYPEDEF 1

	/*============================================================*\
	||                                                            ||
	|| Enthaelt die noetigen Typendeklaraionen fuer Langzahl-     ||
	|| arithmetik, Arithmetik in Fp, Fq, abelschen Zahlkoerpern   ||
	|| und p_adische Koerper.                                     ||
	||                                                            ||
	\*============================================================*/

typedef long polynom_Fp;
typedef long padic;
typedef long Fq;
typedef unsigned long lang;
typedef int boolean;
typedef struct { unsigned long low, hi;	} prod; 
typedef struct { lang st, nd, rd;		} hipos;
typedef struct { int z; int n;			} bruch;
typedef struct { lang z; lang n;		} rational;
typedef struct { lang f1, f2, g1, g2;   } pair;
typedef struct { int  f1, f2, g1, g2;   } spair;

typedef struct { int len; char *string;	} padicchar;
typedef struct { int len; char *string;	} langchar;
typedef struct { int len; char *string;	} Fqchar;

	/*============================================================*\
	||                                                            ||
	|| Definition und Deklaration der Additions- und Multiplika-  ||
	|| tionstabelle der Fp-Arithmetik.                            ||
	|| (schneller, aber nicht kompatibel zur Fq-Arithmetik)       ||
	||                                                            ||
	\*============================================================*/

#define SHORTPRIME_LIMIT 300   /* Grenze fuer Rechnen mit Tabelle  */
#define MIDDLEPRIME_LIMIT 8192     /* Grenze fuer Cech-Logarithmen */
#define LONGPRIME_LIMIT 1000000000 /* Obergrenze fuer Arithmetik */

typedef struct {
	int s_number, *s_primes; /* Anzahl kleiner Pz; Liste der ... */
	int l_number, *l_primes; /* Anzahl grosser Pz; Liste der ... */
	int *s_kind, *l_kind;    /* Art der Pz: -1 : Primzahl,
								1 : Primzahlpotenz, 0 sonst      */
	int ***ADD, ***MUL;	     /* Rechen-Tabellen fuer kleine Pz.  */
	int **LOG, **EXP;      	 /* p-Logarithmen fuer grosse Pz.    */
	} table;
	
table *prime_table;

int **P_ADD, **P_MUL;               /* aktuelle Rechen_Tabelle */
int *P_LOG, *P_EXP;                 /* aktuelle p-Logarithmen  */
int (*S)(int, int), (*P)(int, int); /* aktuelle Funktion       */
int (*min_mod_factor)(int, int);    /* Verallgemeinerung von a/b */
int act_prime;                 /* aktuelle Primzahl       */
int act_prime_code;            /* Aktuelle Kodierung (s.o.) */

	/*============================================================*\
	||                                                            ||
	|| Definition der Typen fuer die p_adische Arithmetik.        ||
	||                                                            ||
	\*============================================================*/

typedef struct {
	boolean init_done;
	char *name;
	long list_num; /* die eigene Koerperkennung */
	long padic_p, p_halbe;
	long padic_e, padic_f;
	long padic_ef;       /* erspart die Rechnung von e*f */
	long padic_ef_plus1; /* e*f+1 */
	long padic_f_plus1;  /* f+1   */
	long padic_2f;       /* 2*f   */
	long padic_2f_1;     /* 2*f+1 */
	long padic_2e_1;     /* 2*e+1 */
	long padic_e2f_1;    /* e*(2*f-1) */
	long padic_2e_12f_1; /* (2*e-1)*(2*f-1) */
	long padic_cd;       /* cd = Rechentiefe */
	long padic_cd_plus1; /* cd+1 */
	polynom_Fp f_poly;   /* unverzweigtes Polynom */
	polynom_Fp e_poly;   /* Eisensteinpolynom */
	boolean f_poly_is_Conway;  /* Polynom automatisch bestimmt */
	long **f_mat;     /* Reduktionsmatrix der unv. Erweiterung */
	long **e_mat;     /* Reduktionsmatrix der verzw. Erweiterung */
	long **mult_mat;  /* e_mat (x) f_mat , d. h. Kroneckerprodukt */
	padic *padicmulvec;
	long ***diag_f_mat; /*Blockdiagonalmatrix zum Invert. unv. Ausdr. */
	long **inv_mat; /*Diese Matrix wird beim Invertieren gebraucht */
	long **inv_mat_inv; /* Die invertierte inv_mat_matrix */
	long *a_inv_vec; /* Zum Invertieren,falls padic_e>1 */
	long *b_inv_vec;  /* " */
	long *c_inv_vec;  /* " */
	long *a_vec;    /* Zum Invertieren, falls padic_f>1 */
	long *b_vec;      /* " */
	long *c_vec;      /* " */
	long *cb0_vec;    /* " */
	long **a_mat; /* Zum Invertieren, falls f>1 und e>1 */
	long **b_mat;     /* " */
	long **c_mat;     /* " */
	long *sum_vec;    /* " */
	padic y_inversp;   /* wird zum schiften bei div gebraucht, falls e>1
		        y_inversp erhaelt den Wert y^(-1)*p */
	padic ypsilon;
        lang *modvec;
        boolean GEMISCHT; /* Flag: TRUE = echt gemischte Erweiterung, d. h.
                                   das Eisensteinpolynom steht in eisen_mat
                                   da es verzweigte Koeffizienten hat */
        long **eisen_mat; /* Das Eisensteinpolynom im gemischten Fall 
                             Der Matrixeintrag eisen_mat[i][j] ist der
                             Koeffizient e_{i,j}*x^j*y^i des Eisensteinpolys */ 
    } padic_field_TYP;

	/*============================================================*\
	||                                        					  ||
	|| Typendefinition fuer die Arithmetik in endlichen Koerpern. ||
	||				     		                                  ||
	\*============================================================*/

typedef struct {
	char *name;  /* Der eigene Name als Matrixkopf: F[p,f,...] */
	long list_num; /* die eigene Koerperkennung */
	boolean init_done_Fq;
	long Fq_p,  Fq_p_halbe;
	long Fq_f,  Fq_f_plus1;
	long Fq_2f, Fq_2f_1;
	long Fq_q, Fq_q_1; /* Anzahl der Elemante, -1 in Fq */
	Fq Fqmulvec;
	long **f_mat_Fq;
	long ***diag_f_mat_Fq;
	long *a_vec_Fq;
	long **inv_mat_Fq; /*Diese Matrix muss beim Invert.invert. werden */
	long **inv_mat_inv_Fq; /* Die invertierte inv_mat_matrix */
	polynom_Fp f_poly_Fq;
	boolean listtyp;           /* Flag ueber Art der Multiplikation */
	boolean f_poly_is_Conway;  /* Polynom automatisch bestimmt */
        Fq   *Fq_list_etoz; /* speichert n -> a  fuer Fq_pr^n=a */
        long *Fq_list_ztoe; /* speichert lex(a) -> n fuer Fq_pr^n=a */
        Fq   Fq_pr;         /* Ein primitives Element */
    } finite_field_TYP;


	/*============================================================*\
	||		  				                                      ||
	|| Typedefinitions for the cyclotomic fields and the Galois-  ||
	|| groups of those (sub-)fields.		                      ||
	||					                                       	  ||
	|| There is no special type-definition for the elements of a  ||
	|| cyclotomic field, and the cyclotomic numbers which are     ||
	|| integers are represented as 'lang'. The others are vectors ||
	|| of integers with the following structure:                  ||
	||                                                            ||
	||     /------------------------------------------\           ||
	||     | Code | field+base| kgv | c_1 | ... | c_n |	          ||
	||     \------------------------------------------/	          ||
	||	    ^identification   ^                                   ||
	||	     number for       coefficients for basis              ||
	||	     field & base     vectors.                            ||
	||                                                            ||
	|| The code number has the following structure:               ||
	||                  code-bits from long-integer code \        ||
	||                                           |                ||
	||                 identifier for cyclotomics |      |        ||
	||    /------------------------------------------------\      ||
	||    |                	length        | 0 s | 1 0 | 0 0|      ||
	||    \------------------------------------------------/      ||
	||    vector length = rank of field       |	                  ||
	||    s = 1: coefficients are integers  --/	                  ||
	||                                                            ||
	\*============================================================*/

	/*============================================================*\
	|| Typedefinition for elements of a Galoisgroup. The elements ||
	|| are given by the image (= exponent) of e_n, and the order  ||
	|| of the automorphism.			         ||
	\*============================================================*/

typedef struct { int image, order; } Galois_element;

	/*============================================================*\
	|| Typedefinition for a abelian Galoisgroup. It contains the  ||
	|| order, the order of the cyclotomic field, and the genera-  ||
	|| tors in normal form.                                       ||
	\*============================================================*/

typedef struct { int group_order, field_order;
				 int gen_num;
				 Galois_element *gen; } Galois_group;

	/*============================================================*\
	|| Typedefinition for the equivalence-classes Ci in [TB].     ||
	|| It contains the number of classes, the classes, sorted in  ||
	|| ascending order of d_i, and for each d_i the standard      ||
	|| primitive roots (indexed).                            	  ||
	\*============================================================*/

typedef struct { int num_class;
				 int **classes; } Ci_class;

	/*============================================================*\
	|| Summand of a base element, contains exponent and	          ||
	|| coefficient of the root.                  	    	      ||
	\*============================================================*/

typedef struct s_summand summand;
struct s_summand { int exponent;
				 lang coeff; 
				 summand *next;};

	/*============================================================*\
	|| Entry of a multiplication table, contains the row, col and ||
	|| coefficient. The definition is recursiv, so the muliplica- ||
	|| tion matrix can be stored as a recursiv list.              ||
	\*============================================================*/
typedef struct s_entry table_entry;
struct s_entry { int col;
				 lang coeff; 
				 table_entry *next;} ;

	/*============================================================*\
	|| Typedefinition for a subfield-base. It contains the order  ||
	|| of a primitive root, the number and the list of base-      ||
	|| vectors, given as linear combinations of the roots, and    ||
	|| the multiplication tables for this base. The coefficients  ||
	|| in the base and in the tables can be cyclotomic numbers.   ||
	|| In the list substitute the linear combinations for those   ||
	|| roots which are not basic elements are given (if possible).||
	|| Here the exponents refers to the basic vectors, not to the ||
	|| root exponents.                                            ||
	|| The i-th entry in the multiplication table points to the   ||
	|| multiplication table for the i-th coefficient of the       ||
	|| product.                                                   ||
	\*============================================================*/

typedef struct { int root_order, rank;
	       prod *order_primes;
				 summand *base, *substitute;
	       table_entry **mult_table; } field_BASE;

typedef struct s_field cyclo_field_TYP;

	/*============================================================*\
	|| Allgemeiner Typ fuer Koerper; dient zur globalen Verwaltung||
	\*============================================================*/

typedef struct {
	int Typ, list_num;
	char *name;
	 cyclo_field_TYP *CF;
	 padic_field_TYP *PF;
	finite_field_TYP *FF;
	} field_TYP;

	/*============================================================*\
	|| Weitere aufbauende Definitionen:                           ||
	|| Definition des matrix_TYP`s,                               ||
	||                                                            ||
	\*============================================================*/
typedef struct {
	lang     **Z, **N;
	int      **SZ, **SN;
	table_entry **ZC;
	} array_TYP;

typedef struct {
   boolean Typ ,
	 Integral ,
		   Symmetric,
		   Diagonal ,
		   Scalar    ,
		   Sparse,
		   Permutation,
		   Monomial;
 	int    Max_Entry ; /* Maximaler Eintrag in short-integer-Matrizen*/
	} flag_TYP;

typedef struct {
	flag_TYP flags;
	int cols, rows, prime;
	lang kgv;
	char *name;      /* Name des Koerpers */
	int field_num;   /* Nummer des Koerpers in allgemeiner Liste */
   	array_TYP array;
	} matrix_TYP;

	/*============================================================*\
	|| Typedefinition for a cyclotomic field. It contains vectors ||
	|| for bases, subfields, supfields and those created field    ||
	|| which are not sub/sup-field. Furthermore it contains the   ||
	|| Galoisgroup, and, if the field is a subfield of a cyclo-   ||
	|| tomic field, it contains a pointer to the cyclotomic field ||
	|| (conductor), and the group which fixes this field.         ||
	|| Often the conductor field is not need arithmetically, so   ||
	|| there is also a pointer to the conductor_base.             ||
	||                                                            ||
	|| The matrix (m)_{ij} contains the base-change coefficient   ||
	|| from the i-th base to the j-th base.                       ||
	||                                                            ||
	|| There is also space for various flags.                     ||
	||                                                            ||
	\*============================================================*/

struct s_field { int list_num, base_num;
				 int supfield_num, subfield_num, parfield_num;
	       int root_order, rank, flags;
				 char *name;
				 field_BASE **bases, *con_base;
				 field_BASE **rel_bases;
				 matrix_TYP **base_change; 
				 matrix_TYP **norms;
				 Ci_class *base_classes;
				 field_TYP **supfields, **subfields, **parfields;
				 field_TYP *conductor;
				 Galois_group *Act_group, *Fix_group;
				 };

	/*============================================================*\
	||                                                            ||
	|| Typendefinition fuer Polynome.                             ||
	|| fac_num enthaelt die Anzahl der gefundenen Faktoren.       ||
	|| Diese werden in factors abgelegt.                          ||
	|| *.hi enthealt jeweils den Exponenten, *.low den entspre-   ||
	|| chenden Koeffizienten.                                     ||
	|| Die Exponenten sind absteigend geordnet.                   ||
	|| 'factorize' ist 'TRUE', falls alle Faktoren irreduzibel    ||
	|| sind, also das Polynom vollstaendig faktorisiert ist.      ||
	||                                                            ||
	\*============================================================*/
typedef struct {
	int fac_num, grad;
	boolean factorize;
	prod **factors;
	} polynom;

typedef struct {
    long  deg;   /* Grad des Polynoms */
    padic *kof;  /* Die Koeffizenten des Polynoms, wobei
                     kof[0] den x^0-Term usw. beinhaltet */
    } padicpoly;

	/* Definition fuer das Arbeiten mit Worten */
typedef struct s_wordchain word_chain;

struct s_wordchain { word_chain *prev, *next;
				 int name, coeff;};
#endif
