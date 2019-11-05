#include "typedef.h"
#include "matrix.h"
#include "datei.h"
#include "getput.h"
#include "bravais.h"
#include "voronoi.h"



int main(int argc, char *argv[])
{
   bravais_TYP *G, *Gtr;
   matrix_TYP *A, *P, *M, *S;
   int i,j;
   int *av, *pv;
   int Fdim;
   int rc, lc, g;

        read_header(argc, argv);
        if(FILEANZ != 3)
        { printf("Usage: Perfect_neighbours 'file1' 'file2' 'file3'\n");
          printf("with: file1 contains the bravais_TYP G\n");
          printf("      file2 contains the bravais_TYP G^{tr}\n");
          printf("      file3 contains a  G-perfect form\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }
        G = get_bravais(FILENAMES[0]);
        Gtr = get_bravais(FILENAMES[1]);
        A = get_mat(FILENAMES[2]);
        Fdim = G->form_no;
        S = trace_bifo(G->form, Gtr->form, Fdim);
        M = all_voronoi_neighbours(A, G, Gtr->form, S);
        if( (av = (int *)malloc(Fdim *sizeof(int))) == NULL)
        {
          printf("malloc of av in main programm failed\n");
          exit(2);
        }
        if( (pv = (int *)malloc(Fdim *sizeof(int))) == NULL)
        {
          printf("malloc of pv in main programm failed\n");
          exit(2);
        }
        form_to_vec(av, A, G->form, Fdim, &i);
        printf("#%d\n", M->rows);
        for(i=0;i<M->rows;i++)
        {
           lc = M->array.SZ[i][Fdim];
           rc = M->array.SZ[i][Fdim+1];
           g = M->array.SZ[i][Fdim+2];
           for(j=0;j<Fdim;j++)
             pv[j] = (lc * av[j] + rc * M->array.SZ[i][j])/g;
           P = vec_to_form(pv, G->form, Fdim);
           put_mat(P, NULL, "G-perfect form", 2);
           free_mat(P);
        }
        printf("Coordinates of initial form:\n");
        for(i=0;i<Fdim;i++)
          printf("%d  ", av[i]);
        printf("\n");
        put_mat(M, NULL, "further informations of other perfect forms", 0);


   exit(0);
}
