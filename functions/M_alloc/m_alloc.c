#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	unsigned long pt;
	long noc;
	} STAT_LIST_TYPE;
int SCOUNT = 0, SFLAG = 0;
int SHELP = 0;
long SMEMORY = 0x10009b50L;
int PERIOD = 1, PCOUNT = 0, P_ACT = 0;
int LIST_SIZE = 0, LIST_USED = 0;
unsigned **MALLOC_LIST = NULL;

int *STAT_LIST_SIZE = NULL;
int *STAT_LIST_USED = NULL;
STAT_LIST_TYPE **STAT_MALLOC_LIST = NULL;

#define M_ALLOC_MAGIC -4711

	/*============================================================*\
	||                                                            ||
	|| Some tools for finding errors in memory allocation:        ||
	||                                                            ||
	\*============================================================*/

	/*============================================================*\
	||                                                            ||
	|| Bilde eine Liste aller Pointer, die allociert wurden.      ||
	|| Die Liste enthaelt den Pointer selbst und einen Zaehler,   ||
	|| wie oft der Pointer allociert und wieder freigegeben wurde.||
	|| Aktive Pointer haben einen positiven Zaehlwert, freigege-  ||
	|| bene einen negativen. Mehrfaches Freigeben fuehrt zum      ||
	|| abbruch.                                                   ||
	||                                                            ||
	\*============================================================*/

void 
pointer_statistics (
    unsigned *p,
    int status  /* 1: pointer wurde neu allociert
				2: pointer wurde freigegeben
				0: Ausgabe der Liste           */
)
{
int i,j;
int flag;
FILE *outfile;
STAT_LIST_TYPE *HELP_MALLOC_LIST;
unsigned long p_cut;

if(status == 0) { /* Gibt die Liste auf eine Datei aus */
	outfile = (FILE *)fopen("Malloc_liste", "w");
	flag = 0;
	fprintf(outfile,"Liste der allocierten und nicht wieder freigegebenen Pointer:\n");
	if(STAT_LIST_USED == NULL) {
		fprintf(outfile, "Keine Pointer mehr aktiv!! \n");
		fclose(outfile);
		return;
		}
	for(j = 0; j < STAT_LIST_USED[-1]; j++) {
		HELP_MALLOC_LIST = STAT_MALLOC_LIST[j];
		for(i = 1; i < STAT_LIST_USED[j]; i++) {
			if(HELP_MALLOC_LIST[i].noc > 0) {
				flag = 1;
				fprintf(outfile, "%lx: %ld\n", HELP_MALLOC_LIST[i].pt,
				   HELP_MALLOC_LIST[i].noc);
				}
			}
		}
	if(flag == 0) {
		fprintf(outfile, "Keine Pointer mehr aktiv!! \n");
		}
	fclose(outfile);
	}

else if(status == 1) { /* Neuen Pointer einfuegen */
	if(STAT_MALLOC_LIST == NULL) { /* Liste anlegen */
		STAT_MALLOC_LIST = (STAT_LIST_TYPE **)malloc(10 * 
											sizeof(STAT_LIST_TYPE *));
		STAT_LIST_USED = (int *)calloc(11, sizeof(int));
		STAT_LIST_SIZE = (int *)calloc(11, sizeof(int));
		STAT_LIST_USED++;
		STAT_LIST_SIZE++;
		STAT_LIST_SIZE[-1] = 10;
		for(i = 0; i < 10; i++) {
			STAT_MALLOC_LIST[i] = (STAT_LIST_TYPE *)calloc(100,
									sizeof(STAT_LIST_TYPE));
			STAT_LIST_SIZE[i] = 100;
			STAT_LIST_USED[i] = 1;
			}
		}
	p_cut = ((unsigned long) p) >> 14; /* Hoffentlich realistische Groesse */
	i = 0;
	while((i < STAT_LIST_USED[-1]) && 
		  (p_cut != STAT_MALLOC_LIST[i][0].pt)) i++;
	if(i < STAT_LIST_USED[-1]){ /* Pointerbereich schon mal allociert */
		HELP_MALLOC_LIST = STAT_MALLOC_LIST[i];
		j = 1;
		while((j < STAT_LIST_USED[i]) && 
		  	(p != (unsigned *)HELP_MALLOC_LIST[j].pt)) j++;
		if(j < STAT_LIST_USED[i]) { /* Pointer schon mal allociert */
			if(HELP_MALLOC_LIST[j].noc > 0) {
				fprintf(stderr,"Error in pointer_statistics.\n");
				}
			HELP_MALLOC_LIST[j].noc = 1 - HELP_MALLOC_LIST[j].noc;
			}
		else {
			if(STAT_LIST_USED[i] == STAT_LIST_SIZE[i]) {
				HELP_MALLOC_LIST = (STAT_LIST_TYPE *)realloc(
					 HELP_MALLOC_LIST, (STAT_LIST_SIZE[i]+100)
						 *sizeof(STAT_LIST_TYPE ));
				STAT_LIST_SIZE[i] += 100;
				STAT_MALLOC_LIST[i] = HELP_MALLOC_LIST;
				}
			HELP_MALLOC_LIST[STAT_LIST_USED[i]  ].pt = (unsigned long)p;
			HELP_MALLOC_LIST[STAT_LIST_USED[i]++].noc = 1;
			}
		}
	else {
		if(STAT_LIST_USED[-1] == STAT_LIST_SIZE[-1]) {
			STAT_LIST_USED--;
			STAT_LIST_USED = (int *)realloc(STAT_LIST_USED,
							       (STAT_LIST_SIZE[-1]+11)*sizeof(int));
			STAT_LIST_USED++;
			STAT_LIST_SIZE--;
			STAT_LIST_SIZE = (int *)realloc(STAT_LIST_SIZE,
							       (STAT_LIST_USED[-1]+11)*sizeof(int));
			STAT_LIST_SIZE++;
			STAT_LIST_SIZE[-1] += 10;
			STAT_MALLOC_LIST = (STAT_LIST_TYPE **)realloc(
							   STAT_MALLOC_LIST, (STAT_LIST_SIZE[-1])
							   *sizeof(STAT_LIST_TYPE *));
			for(j = STAT_LIST_USED[-1]; j < STAT_LIST_SIZE[-1]; j++) {
				STAT_MALLOC_LIST[j] = (STAT_LIST_TYPE *)calloc
						(100,sizeof(STAT_LIST_TYPE));
				STAT_LIST_SIZE[j] = 100;
				STAT_LIST_USED[j] = 1;
				}
			}
		STAT_MALLOC_LIST[i][0                  ].pt = p_cut;
		STAT_MALLOC_LIST[i][STAT_LIST_USED[i]  ].pt = (unsigned long)p;
		STAT_MALLOC_LIST[i][STAT_LIST_USED[i]++].noc = 1;
		STAT_LIST_USED[-1]++;
		}
	}
else if(status == 2) {
	i = 0;
	p_cut = ((unsigned long)p)>> 14;
	while((i<STAT_LIST_USED[-1])&&(p_cut!=STAT_MALLOC_LIST[i][0].pt))
		i++;
	if(i < STAT_LIST_USED[-1]) { /*Pointerbereich schon mal allociert */
		j = 1;
		HELP_MALLOC_LIST = STAT_MALLOC_LIST[i];
		while((j < STAT_LIST_USED[i]) && 
			(p != (unsigned *)HELP_MALLOC_LIST[j].pt))
			j++;
		if(j < STAT_LIST_USED[i]) { /*Pointer schon mal allociert */
			if(HELP_MALLOC_LIST[j].noc > 0) { /* Alles okay */
				HELP_MALLOC_LIST[j].noc *= (-1);
				}
			else {
				fprintf(stderr, 
					"Fehler: Doppeltes Freigeben eines Pointers.\n");
				}
			}
		else {
			fprintf(stderr, 
				"Fehler: Freigeben von nicht allociertem Pointer.\n");
			}
		}
	else {
		fprintf(stderr, 
			"Fehler: Freigeben eines nicht allocierten Pointers.\n");
		}
	}
}

	/*============================================================*\
	||                                                            ||
	|| First version for allocation-error-diagnostics:            ||
	|| Get 4 more bytes before and after the array and write      ||
	|| certain values in it.                                      ||
	|| Check the values when freeing it.                          ||
	||                                                            ||
	\*============================================================*/
int *m_alloc_d1(int size_t)
{
int *p;
int newsize;

newsize = (size_t-1) / 4 + 1;

if( (p = (int*)malloc(4*newsize+32)) == NULL) {
	fprintf(stderr,"Fehler in malloc \n");
	exit(2);
	}
if(SFLAG) pointer_statistics(p, 1);
p[3] = newsize;
p[1] = p[2] = p[0] = M_ALLOC_MAGIC;
p += 4;
p[newsize] = p[newsize+1] = p[newsize+2] = p[newsize+3] = M_ALLOC_MAGIC;
if(p == (int *)SMEMORY) {
	SHELP++;
	}
SCOUNT++;
return(p);
}

int *c_alloc_d1(int size_t, int size_n)
{
int *p;
int newsize;

newsize = (size_t*size_n -1) / 4 + 1;
if( (p = (int*)malloc(4*newsize+32)) == NULL) {
	fprintf(stderr,"Fehler in calloc \n");
	exit(2);
	}
if(SFLAG) pointer_statistics(p, 1);
p[3] = newsize;
p[1] = p[2] = p[0] = M_ALLOC_MAGIC;
p += 4;
p[newsize] = p[newsize+1] = p[newsize+2] = p[newsize+3] = M_ALLOC_MAGIC;
while((--newsize) >= 0) {
	p[newsize] = 0;
	}
if(p == (int *)SMEMORY) {
	SHELP++;
	}
SCOUNT++;
return(p);
}

int *re_alloc_d1(int *old_p, int size_t)
{
int *p;
int oldsize, sizestore;

if(old_p == NULL)
  return(m_alloc_d1(size_t));
 
old_p -= 4;
if((old_p[1]!=M_ALLOC_MAGIC)||(old_p[2]!=M_ALLOC_MAGIC)||(old_p[0]!=M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in re_alloc-control sequence\n");
	exit(2);
	}
oldsize = old_p[3]+4;
if((old_p[oldsize  ] != M_ALLOC_MAGIC) || (old_p[oldsize+1] != M_ALLOC_MAGIC) ||
   (old_p[oldsize+2] != M_ALLOC_MAGIC) || (old_p[oldsize+3] != M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in re_alloc-control sequence\n");
	exit(2);
	}

oldsize -= 4;
sizestore = (size_t-1) / 4 + 1;
if(SFLAG) pointer_statistics(old_p,2);
if( (p = (int*)realloc(old_p, sizestore*4 +32)) == NULL) {
	fprintf(stderr,"Fehler in realloc \n");
	exit(2);
	}
if(SFLAG) pointer_statistics(p,1);
p[3] = sizestore;
p[1] = p[2] = p[0] = M_ALLOC_MAGIC;
p += 4;
if(p == (int *)SMEMORY) {
	SHELP++;
	}
if(sizestore > oldsize) {
	memset(p + oldsize, 0, (sizestore-oldsize) * sizeof(int));
	}
p[sizestore] = p[sizestore+1] = p[sizestore+2] = p[sizestore+3] = M_ALLOC_MAGIC;
return(p);
}

void 
fr_ee_d1 (int *p)
{
int oldsize;

if(p == (int *)SMEMORY) {
	SHELP++;
	}
p -= 4;
oldsize = p[3] + 4;
if((p[1] != M_ALLOC_MAGIC) || (p[2] != M_ALLOC_MAGIC) || (p[0] != M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in fr_ee-control sequence\n");
	exit(2);
	}
if( (p[oldsize  ] != M_ALLOC_MAGIC) || (p[oldsize+1] != M_ALLOC_MAGIC) || 
	(p[oldsize+2] != M_ALLOC_MAGIC) || (p[oldsize+3] != M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in fr_ee-control sequence\n");
	exit(2);
	}
SCOUNT--;
if(SFLAG) pointer_statistics(p, 2);
free(p);
}


	/*============================================================*\
	||                                                            ||
	|| Second version for allocation-error-diagnostics:           ||
	|| Get 4 more bytes before and after the array and write      ||
	|| certain values in it.                                      ||
	|| Create a list of all allocated pointers and check the      ||
	|| values periodically.                                       ||
	||                                                            ||
	\*============================================================*/

void 
add_pointer (unsigned *p)

{
int i, oldsize;

if(MALLOC_LIST == NULL) {
	MALLOC_LIST = (unsigned **)calloc(100,sizeof(unsigned *));
	LIST_SIZE = 100;
	LIST_USED = 0;
	}
if(LIST_SIZE == LIST_USED) {
	MALLOC_LIST = (unsigned **)realloc(MALLOC_LIST,
					(100+LIST_SIZE)* sizeof(unsigned *));
	memset(MALLOC_LIST+LIST_SIZE,0,100*sizeof(unsigned ));
	LIST_SIZE += 100;
	}
while(MALLOC_LIST[P_ACT] != NULL) P_ACT = (P_ACT + 1)%LIST_SIZE;
MALLOC_LIST[P_ACT] = p;
LIST_USED++;
PCOUNT ++;
if(PCOUNT == PERIOD) {
	PCOUNT = 0;
	for(i = 0; i < LIST_SIZE; i++) {
		if(MALLOC_LIST[i] != NULL) {
			p = MALLOC_LIST[i] - 4;
			oldsize = p[3] + 4;
			if((p[1] != M_ALLOC_MAGIC) || (p[2] != M_ALLOC_MAGIC) || (p[0] != M_ALLOC_MAGIC)) {
				fprintf(stderr,"Error in alloc-control sequence\n");
				exit(2);
				}
			if( (p[oldsize  ] != M_ALLOC_MAGIC) || (p[oldsize+1] != M_ALLOC_MAGIC) || 
				(p[oldsize+2] != M_ALLOC_MAGIC) || (p[oldsize+3] != M_ALLOC_MAGIC)) {
				fprintf(stderr,"Error in alloc-control sequence\n");
				exit(2);
				}
			}
		}
	}
}

void 
delete_pointer (unsigned *p)

{
int i, oldsize;

i = 0;
while((i < LIST_SIZE) && MALLOC_LIST[i] != p) i++;
if(i == LIST_SIZE) {
	fprintf(stderr,"Error in alloc-sequence:\n");
	fprintf(stderr,"Free-call on non-allocated pointer\n");
	exit(2);
	}
MALLOC_LIST[i] = NULL;
LIST_USED--;
PCOUNT ++;
p -= 4;
if((p[1]!=M_ALLOC_MAGIC)||(p[2]!=M_ALLOC_MAGIC)||(p[0]!=M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in alloc-control sequence\n");
	exit(2);
	}
oldsize = p[3]+4;
if((p[oldsize  ] != M_ALLOC_MAGIC) || (p[oldsize+1] != M_ALLOC_MAGIC) ||
   (p[oldsize+2] != M_ALLOC_MAGIC) || (p[oldsize+3] != M_ALLOC_MAGIC)) {
	fprintf(stderr,"Error in re_alloc-control sequence\n");
	exit(2);
	}
if(PCOUNT == PERIOD) {
	PCOUNT = 0;
	for(i = 0; i < LIST_SIZE; i++) {
		if(MALLOC_LIST[i] != NULL) {
			p = MALLOC_LIST[i] - 4;
			oldsize = p[3] + 4;
			if((p[1] != M_ALLOC_MAGIC) || (p[2] != M_ALLOC_MAGIC) || (p[0] != M_ALLOC_MAGIC)) {
				fprintf(stderr,"Error in alloc-control sequence\n");
				exit(2);
				}
			if( (p[oldsize  ] != M_ALLOC_MAGIC) || (p[oldsize+1] != M_ALLOC_MAGIC) || 
				(p[oldsize+2] != M_ALLOC_MAGIC) || (p[oldsize+3] != M_ALLOC_MAGIC)) {
				fprintf(stderr,"Error in alloc-control sequence\n");
				exit(2);
				}
			}
		}
	}
}


int *m_alloc_d2(int size_t)
{
int *p;
int newsize;

newsize = size_t / 4 + 1;

if( (p = (int*)malloc(4*newsize+32)) == NULL) {
	fprintf(stderr,"Fehler in malloc \n");
	exit(2);
	}
if(SFLAG) pointer_statistics(p, 1);
p[3] = newsize;
p[1] = p[2] = p[0] = M_ALLOC_MAGIC;
p += 4;
p[newsize] = p[newsize+1] = p[newsize+2] = p[newsize+3] = M_ALLOC_MAGIC;
if(p == (int *)SMEMORY) {
	SHELP++;
	}
SCOUNT++;

add_pointer(p);

return(p);
}

int *c_alloc_d2(int size_t, int size_n)
{
int *p;
int newsize;
if( (p = (int*)m_alloc_d2(size_t*size_n)) == NULL) {
	fprintf(stderr,"Fehler in calloc \n");
	exit(2);
	}
memset(p,0,size_t*size_n);
return(p);
}

int *re_alloc_d2(int *old_p, int size_t)
{
int *p;
int sizestore;

if(old_p == NULL)
  return(m_alloc_d2(size_t));

delete_pointer(old_p);

old_p -= 4;

sizestore = size_t / 4 + 1;
if(SFLAG) pointer_statistics(old_p, 2);
if( (p = (int*)realloc(old_p, sizestore*4 +32)) == NULL) {
	fprintf(stderr,"Fehler in realloc \n");
	exit(2);
	}
if(SFLAG) pointer_statistics(p, 1);
p[3] = sizestore;
p[1] = p[2] = p[0] = M_ALLOC_MAGIC;
p += 4;
p[sizestore] = p[sizestore+1] = p[sizestore+2] = p[sizestore+3] = M_ALLOC_MAGIC;
if(p == (int *)SMEMORY) {
	SHELP++;
	}
add_pointer(p);

return(p);
}

void 
fr_ee_d2 (int *p)
{
int oldsize;

delete_pointer(p);
if(p == (int *)SMEMORY) SHELP++;
p -= 4;
SCOUNT--;
if(SFLAG) pointer_statistics(p, 2);
free(p);
}

