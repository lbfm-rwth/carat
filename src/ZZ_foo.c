#include"typedef.h"
#include"tools.h"
#include"matrix.h"
#include"bravais.h"
#include"ZZ.h"

void scan_argv (argc, argv, filename)
int argc;
char *argv[], **filename;
{
if ( argc != 2 )
	*filename = NULL;
else
	*filename = argv[1];
}

int ABBRUCH;
/*{{{}}}*/
/*{{{  ZZ_usage*/
void ZZ_usage(progname)
char *progname;
{
  fprintf(stderr,"Usage: %s -bghl <#level> n <#number> pqrstu <file>\n\n",progname);
  fprintf(stderr,"-b  : Print only the matrices of change of base and their inverse.\n");
  fprintf(stderr,"-g  : Do not compute elementary divisors of the gram matrix.\n");
  fprintf(stderr,"-h  : Print this help\n");
  fprintf(stderr,"-l #: Stop after reaching level #level (default #=%d).\n",LEVEL);
  fprintf(stderr,"-n #: Stop after computation of #number \"Zentrierungen\" (default #=%d).\n",NUMBER);
  fprintf(stderr,"-p<d0>/<d1>/<d2> ... : treat the lattice as a direct sum of <d0> sublattices\n");
  fprintf(stderr,"      of dimensions <d1>, <d2> etc. (0 <= d0 <= 6) and compute only those\n");
  fprintf(stderr,"      centerings that have surjective projections on them.\n");
  fprintf(stderr,"-q  : Quiet mode. Suppress any messages to stdin/stdout.\n");
  fprintf(stderr,"-r  : With ZZ_lll-reduction.\n");
  fprintf(stderr,"-s  : Print less information.\n");
  fprintf(stderr,"-t  : Create the data-file \"ZZ.tmp\".\n");
  fprintf(stderr,"-u  : Do not compute elementary divisors of the change of base\n\n");
}

/*}}}  */


void foo()
{
/* printf("Hello world!\n"); */
}

void main (argc, argv)

int argc;
char *argv[];
{  
matrix_TYP *Gram, **help2;
ZZ_data_t data;
ZZ_tree_t tree;

int i, j, k;

char *file_name;

  scan_argv(argc, argv, &file_name);
  Gram= ZZ_fget_data (&data,&tree,file_name);
  if ( constituents == 1)
  {
          for ( i = 0; i < data.p_consts.k; i++)
          {
                  help2= ZZ_irr_const( data.DELTA, data.r,
                                          data.p_consts.p[i],
                                          &data.p_consts.s[i]);
            data.n[i] = (int *)malloc(data.p_consts.s[i]*sizeof(int));
                  data.p_consts.Delta[i]=(matrix_TYP ***)malloc(data.p_consts.s[i]*sizeof(matrix_TYP **));
            for ( j = 0 ; j < data.p_consts.s[i]; j ++ )
            {
              data.p_consts.Delta[i][j]= (matrix_TYP **)malloc(data.r*sizeof(matrix_TYP *));
              for ( k = 0; k < data.r; k ++ ) {
                data.p_consts.Delta[i][j][k] = help2[j*data.r+k];
                data.p_consts.Delta[i][j][k]->prime = data.p_consts.p[i];
              }
              data.n[i][j] = data.p_consts.Delta[i][j][0]->rows;
            }
            free( help2 );
          }
          ZZ_test_konst (&data);
                /*{{{  */
                        /*{{{  */
                        for ( i = 0; i < data.p_consts.k; i++) {
                                /*{{{  */
                                /*------------------------------------------------------------*\
                                | initialize Endomorphisms |
                                \*------------------------------------------------------------*/
                                data.EnCo[i]  = (ZZ_prod_t *)malloc(data.p_consts.s[i]*sizeof(ZZ_prod_t));
                                data.Endo[i] =(matrix_TYP ***)malloc(data.p_consts.s[i]*sizeof(matrix_TYP **));
                                  }

                                data.epi_base = NULL;
                                data.epi = init_mat(data.N,data.N,"ik");
                                tree.root->k_vec = (int **)malloc(data.p_consts.k*sizeof(int *));
                                data.VK =(int **)malloc(data.p_consts.k*sizeof(int *));
                                for(i = 0; i < data.p_consts.k; i++) {
                                  tree.root->k_vec[i] = (int *)calloc(data.p_consts.s[i],sizeof(int));
                                  data.VK[i] = (int *)calloc(data.p_consts.s[i]+1,sizeof(int));
                                  data.VK[i]++;
                                  }
                                ZZ_make_endo (&data);
                                /*}}}  */
                        /*}}}  */
                /*}}}  */
  }
  if( verbose == TRUE )
  {
        /*{{{  */
          for ( i = 0; i < data.p_consts.k; i++)
          {
            fprintf(stderr,"Primzahl: %d\n",data.p_consts.p[i]);
            for ( j = 0 ; j < data.p_consts.s[i]; j ++ )
            {
              fprintf(stderr,"Konstituent %d:\n",j);
              for ( k = 0; k < data.r; k ++ )
              {
                fprintf(stderr,"Erzeuger %d:\n",k);
                fput_mat ( stderr, data.p_consts.Delta[i][j][k], "Konstituent", 0);
              }
            }
          }
        /*}}}  */
        }
  ZZ_intern( Gram, &data, &tree );
  ZZ_fput_data (&data, &tree,ABBRUCH);
  ZZ_free_data(&data);
  free(Gram);
}

