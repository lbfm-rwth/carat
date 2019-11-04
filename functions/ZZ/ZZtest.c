#ifdef TEST_PRIVATE
#include "ZZ_P.h"
#include "ZZ_irr_const_P.h"
#include "ZZ_cen_fun_P.h"
#else
#include "ZZ.h"
#endif

/*{{{}}} */
/*{{{  ZZ_usage */
void 
ZZ_usage (char *progname)
{
  fprintf (stderr, "Usage: %s -bghl <#level> n <#number> pqrstu <file>\n\n", progname);
  fprintf (stderr, "-b  : Print only the matrices of change of base and their inverse.\n");
  fprintf (stderr, "-g  : Do not compute elementary divisors of the gram matrix.\n");
  fprintf (stderr, "-h  : Print this help\n");
  fprintf (stderr, "-l #: Stop after reaching level #level (default #=%d).\n", LEVEL);
  fprintf (stderr, "-n #: Stop after computation of #number \"Zentrierungen\" (default #=%d).\n", NUMBER);
  fprintf (stderr, "-p<d0>/<d1>/<d2> ... : treat the lattice as a direct sum of <d0> sublattices\n");
  fprintf (stderr, "      of dimensions <d1>, <d2> etc. (0 <= d0 <= 6) and compute only those\n");
  fprintf (stderr, "      centerings that have surjective projections on them.\n");
  fprintf (stderr, "-q  : Quiet mode. Suppress any messages to stdin/stdout.\n");
  fprintf (stderr, "-r  : With ZZ_lll-reduction.\n");
  fprintf (stderr, "-s  : Print less information.\n");
  fprintf (stderr, "-t  : Create the data-file \"ZZ.tmp\".\n");
  fprintf (stderr, "-u  : Do not compute elementary divisors of the change of base\n\n");
}

/*}}}  */


void 
foo (void)
{
/* printf("Hello world!\n"); */
}

#ifndef TEST_PRIVATE
/*{{{  main */
void 
main (int argc, char *argv[])
{
  bravais_TYP *group;
  matrix_TYP *gram;
  char *file_name = NULL;
  char *prog_name;
  int retval;
  char options[256];
  int is_option;
  int projections[7];
  char *help;
  int num_proj;

  /* skipped the function atexit because some compiler won't have it
  tilman 28/07/97:
  atexit (foo); */
  prog_name = argv[0];
  options[0] = '\0';
  is_option = 0;
  argv++;
  argc--;
  while (argc >= 2)
    {
      if (**argv == '-')
	{
	  (*argv)++;
	  is_option = 1;
	}
      if (is_option)
	{
	  strcat (options, *argv);
	  if ((help = strchr (*argv, 'p')) != NULL)
	    {
	      projections[0] = atoi (help + 1);
	      if (projections[0] == 0)
		{		/* is in next argv */
		  argv++;
		  argc--;
		  projections[0] = atoi (argv[0]);
		}
	      if (projections[0] == 0)
		{
		  fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
		  ZZ_usage (prog_name);
		  exit (31);
		}
	      else if (projections[0] > 6)
		{
		  fprintf (stderr, "Maximal dimensionis 6\n");
		}
	      else
		{
		  printf ("%d\n", projections[0]);
		  help = argv[0];
		  for (num_proj = 1; num_proj <= projections[0]; num_proj++)
		    {
		      if ((help = strchr (help, '/')) != NULL)
			{
			  help++;
			  projections[num_proj] = atoi (help);
			  if (projections[num_proj] == 0)
			    {
			      fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
			      ZZ_usage (prog_name);
			      exit (31);
			    }
			}
		      else
			{
			  fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
			  ZZ_usage (prog_name);
			  exit (31);
			}
		    }
		}
	    }
	}
      argv++;
      argc--;
    }
  file_name = *argv;
  if (file_name == NULL)
    {
      fprintf (stderr, "No filename specified.\n");
      ZZ_usage (prog_name);
      exit (31);
    }

  group = get_bravais (file_name);

  gram = group->form[0];

  ZZ (group, gram, group->divisors, options, projections[0],
      projections[1],
      projections[2],
      projections[3],
      projections[4],
      projections[5],
      projections[6]);

  fprintf (stderr, "num_zentr: %d\n", group->zentr_no);
  cleanup_prime ();
  free_bravais (group);

#ifdef MEM_DEBUG
  debug_meminfo ();
#endif

  exit (retval);

}
/*}}}  */
#else
/*{{{  main      */
void 
main (int argc, char *argv[])
{
  matrix_TYP *Gram, **help2;
  ZZ_data_t data;
  ZZ_tree_t tree;

  int i, j, k;

  char *file_name;

  scan_arg (argc, argv, &file_name);
  Gram = ZZ_fget_data (&data, &tree, file_name);
  if (constituents == 1)
    {
      for (i = 0; i < data.p_consts.k; i++)
	{
	  help2 = ZZ_irr_const (data.DELTA, data.r,
				data.p_consts.p[i],
				&data.p_consts.s[i]);
	  data.n[i] = (int *) malloc (data.p_consts.s[i] * sizeof (int), "main:data.n[i]");
	  data.p_consts.Delta[i] = (matrix_TYP ***) malloc (data.p_consts.s[i] * sizeof (matrix_TYP **), "main:data.p_consts.Delta[i]");
	  for (j = 0; j < data.p_consts.s[i]; j++)
	    {
	      data.p_consts.Delta[i][j] = (matrix_TYP **) malloc (data.r * sizeof (matrix_TYP *), "main:data.p_consts.Delta[i][j]");
	      for (k = 0; k < data.r; k++)
		{
		  data.p_consts.Delta[i][j][k] = help2[j * data.r + k];
		  data.p_consts.Delta[i][j][k]->prime = data.p_consts.p[i];
		}
	      data.n[i][j] = data.p_consts.Delta[i][j][0]->rows;
	    }
	  free (help2);
	}
      ZZ_test_konst (&data);
      /*{{{  */
      /*{{{  */
      for (i = 0; i < data.p_consts.k; i++)
	{
	  /*{{{  */
				/*------------------------------------------------------------*\
                                | initialize Endomorphisms |
                                \*------------------------------------------------------------*/
	  data.EnCo[i] = (ZZ_prod_t *) malloc (data.p_consts.s[i] * sizeof (ZZ_prod_t), "main:data.EnCo[i]");
	  data.Endo[i] = (matrix_TYP ***) malloc (data.p_consts.s[i] * sizeof (matrix_TYP **), "main:data.Endo[i]");
	}

      data.epi_base = NULL;
      data.epi = init_mat (data.N, data.N, "ik");
      tree.root->k_vec = (int **) malloc (data.p_consts.k * sizeof (int *), "main:tree.root->k_vec");
      data.VK = (int **) malloc (data.p_consts.k * sizeof (int *), "main:data.VK");
      for (i = 0; i < data.p_consts.k; i++)
	{
	  tree.root->k_vec[i] = (int *) calloc (data.p_consts.s[i], sizeof (int), "main:tree.root->k_vec[i]");
	  data.VK[i] = (int *) calloc (data.p_consts.s[i] + 1, sizeof (int), "main:data.VK[i]");
	  data.VK[i]++;
	}
      ZZ_make_endo (&data);
      /*}}}  */
      /*}}}  */
      /*}}}  */
    }
  if (verbose == TRUE)
    {
      /*{{{  */
      for (i = 0; i < data.p_consts.k; i++)
	{
	  fprintf (stderr, "Primzahl: %d\n", data.p_consts.p[i]);
	  for (j = 0; j < data.p_consts.s[i]; j++)
	    {
	      fprintf (stderr, "Konstituent %d:\n", j);
	      for (k = 0; k < data.r; k++)
		{
		  fprintf (stderr, "Erzeuger %d:\n", k);
		  fput_mat (stderr, data.p_consts.Delta[i][j][k], "Konstituent", 0);
		}
	    }
	}
      /*}}}  */
    }
  ZZ_intern (Gram, &data, &tree);
  ZZ_fput_data (&data, &tree);
  ZZ_free_data (&data);
  free (Gram);
}
/*}}}  */
#endif
