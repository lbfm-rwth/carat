#include "typedef.h"
#include "matrix.h"
#include "bravais.h"
#include "orbit.h"
#include "getput.h"
#include "sort.h"
#include "datei.h"
int SFLAG;
int INFO_LEVEL;

int main (int argc, char *argv[])
{

	matrix_TYP **M,
                    *Mat,
                   **erg,
                   **representatives,
                   **tmp;

        bravais_TYP *G,
                    *Stab,
                    *UNITY;   /* the unity group */

        int length,
            l,
            i,
            j,
            no_of_orbits=1,
            Manz;

        int *option;

        extern char **FILENAMES;
        extern int FILEANZ;

        read_header(argc, argv);

	if (is_option('h')){
		INFO_LEVEL = optionnumber('h');
	}
	if (INFO_LEVEL == 8) SFLAG = 1;

        if(FILEANZ != 2)
        {
          printf("Usage: %s 'file1' 'file2' [-i] [-r] [-l] [-k] [-t] [-L=n] [-S=n] [-p] [-u] [-g] [-R]\n",argv[0]);
          printf("\n");
          printf("file1:  matrix_TYP, contains a matrix X whose orbit is to be calculated \n");
          printf("file2:  bravais_TYP, contains generators of a group G\n");
          printf("\n");
          printf("Calulates the orbit of the matrix X in file1 under the group G in file2,\n");
          printf("where the action is specified by the options. Default option is action by\n");
          printf("left multiplication.\n");
          printf("\n");
          printf("Options:\n");
          printf("-i     : Use the generators given in file2 and their\n");
          printf("         inverses to calculate the orbit.\n");
          printf("-r     : Operate from the right.\n");
          printf("-l     : Operate from the left (default).\n");
          printf("-k     : Operate via conjugation, ie.  x -> g x g^-1 \n");
          printf("-L=n   : Calculate at most n elements of the orbit.\n");
          printf("         0 means infinity.\n");
          printf("-S=n   : If given as -S or -S=0 a generating set for\n");
          printf("         the stabilizer is calculated. If given as\n");
          printf("         -S=n at most n matrices of the stabilizer\n");
          printf("         are calculated.\n");
          printf("         S=-1 means ONLY the stabilizer is calculated.\n");
          printf("-p     : Operate on pairs of the form {M,-M}.\n");
          printf("-u     : Operate on the set of rows of the matrix given\n");
          printf("         in file1.\n");
          printf("-f     : Operates on quadratic forms via x -> g^-tr x g^-1\n");
          printf("-g     : Operate on sublattices of Z^n spanned by the columns\n");

          printf("          (rows) of the matrices gX (Xg) with g in G. Brakets \n");
          printf("         apply if given with the -r option.\n");
          printf("-R     : Give representatives of the G-orbits at the end of\n");
          printf("         the output.\n");
          printf("WARNING: If the orbit is infinite use option -L!\n");
          printf("\n");
          printf("Cf. Order, Is_finite.\n");
          /* printf("Cf. Orbit_representatives, Order, Is_finite.\n"); */
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }

        M = mget_mat(FILENAMES[0],&Manz);
	Mat = M[0];
        G = get_bravais(FILENAMES[1]);
        option = make_orbit_options();
        Stab = (bravais_TYP *)malloc(sizeof(bravais_TYP));
        Stab->gen_no = Stab->form_no = Stab->zentr_no = 0;
        Stab->order = Stab->normal_no = Stab->cen_no = 0;
        Stab->divisors[0]=1;
        UNITY = init_bravais(G->dim);
        UNITY->gen = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP *));
	UNITY->gen[0] = init_mat(G->dim,G->dim,"1");
        UNITY->gen_no = 1;
        representatives = (matrix_TYP **) malloc(Manz * sizeof(matrix_TYP *));

        if (is_option('S') && Manz > 1){
           fprintf(stderr,"which matrix to you want to calculate the\n");
           fprintf(stderr,"stabilizer of?");
           exit(3);
        }

        erg = orbit_alg(Mat, G, Stab, option, &length);
        representatives[0] = erg[0];

        for (i=1;i<Manz;i++){

            mat_quicksort(erg,0,length-1,mat_comp);

            /* standartize M[i] (in a funny way) */
            tmp = orbit_alg(M[i], UNITY, Stab, option, &l);
            if (l != 1){
               fprintf(stderr,"ERROR in orbit\n");
               exit(3);
            }
            free_mat(M[i]);
            M[i] = tmp[0];
            free(tmp);

            if (mat_search(M[i],erg,length,mat_comp) == -1){
               tmp = orbit_alg(M[i],G, Stab, option, &l);
               erg = (matrix_TYP **) realloc(erg,(length+l)*sizeof(matrix_TYP));
               for (j=0;j<l;j++){
                  erg[length+j] = tmp[j];
               }
               length += l;
               representatives[no_of_orbits] = tmp[0];
               no_of_orbits++;
               free(tmp);
            }
        }

        if (!is_option('S') || optionnumber('S') >= 0){
            printf("#%d\n", length);
            for(i=0;i<length;i++)
               put_mat(erg[i], NULL, "", 2);
        }
        if(is_option('S') == TRUE)
           put_bravais(Stab, NULL, "Stabilizer of the operation");

        if (is_option('R')){
            printf("===== Number of orbits: \n#%d\n",no_of_orbits);
            for (i=0;i<no_of_orbits;i++){
               put_mat(representatives[i],0,"representative of orbit",2);
            }
        }

	free_bravais(UNITY);
	free_bravais(Stab);
	free_bravais(G);
        for(i=0;i<length;i++)
           free_mat(erg[i]);

        for (i=0;i<Manz;i++)
           free_mat(M[i]);
        free(M);
	free(erg);
        free(representatives);
        if (option != NULL) free(option);
	if (INFO_LEVEL == 8) pointer_statistics(0,0);

   exit(0);
}
