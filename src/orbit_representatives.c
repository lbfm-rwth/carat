#error THIS FILE IS NOT USED

#include"typedef.h"
#include"orbit.h"
#include"getput.h"
#include"sort.h"

main (int argc, char *argv[])
{

	matrix_TYP **Mat, *erg;
        bravais_TYP *G;
        int i, Manz, length;
        int *option;

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
           printf("Usage:  %s 'file1' 'file2' [-i] [-r] [-l] [-k] [-t] [-L=n] [-S=n] [-p] [-u]\n",argv[0]);
           printf("\n");
           printf("file1:  matrix_TYP, contains matrices X_i forming a G-set (union of G-orbits)\n");
           printf("file2:  bravais_TYP, contains generators of a group G\n");
           printf("\n");
           printf("Calulates the orbit representatives of the G-orbits on {X_1,...},\n");
           printf("where the action is specified by the options. Default option is action by\n");
           printf("left multiplication.\n");
           printf("\n");
           printf("Options:\n");
           printf("-i     : Use the generators given in file2 and their\n");
           printf("         inverses to calculate the orbit.\n");
           printf("-r     : Operate from the right.\n");
           printf("-l     : Operate from the left (default).\n");
           printf("-k     : Operate via conjugation, ie.  x -> g x g^-1 \n");
           printf("-L=n   : Calculate at most n elements of the Orbit.\n");
           printf("         0 means infinity.\n");
           printf("-S=n   : If given as -S or -S=0 a generating set for\n");
           printf("         the stabilizer is calculated. If given as\n");
           printf("         -S=n at most n matrices of the stabilizer\n");
           printf("         are calculated.\n");
           printf("-p     : Operate on pairs of the form {M,-M}.\n");
           printf("-u     : Operate on the set of rows of the matrix given\n");
           printf("         in file1.\n");
           printf("-g     : Operate on sublattices of Z^n spanned by the columns of the \n");
           printf("         matrices gX with g in G.\n");
           printf("-f     : Operates on quadratic forms via x -> g^-tr x g^-1\n");
           printf("\n");
           printf("Cf. Orbit.\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = mget_mat (FILENAMES[0], &Manz);
        G = get_bravais(FILENAMES[1]);
        option = make_orbit_options();
        mat_quicksort(Mat, 0, Manz-1, mat_comp);
        erg = orbit_representatives(Mat, Manz, G, option, &length, 1);
        printf("#%d\n", length);
        for(i=0;i<length;i++)
           put_mat(Mat[erg->array.SZ[0][i]], NULL, "", 2);
        printf(" length of i-th orbit is\n");
        for(i=0;i<length;i++)
          printf("%d,  ", erg->array.SZ[1][i]);
        printf("\n");

   exit(0);
}
