#include "typedef.h"
#include "contrib.h"
#include "matrix.h"
#include "getput.h"
#include "idem.h"
#include "sort.h"
#include "orbit.h"
#include "bravais.h"
#include "datei.h"

static int *trace_list_B;
static int *order_list_B;

static int *traces_of_inv_A;
static int **traces_of_prod_A;

static int *image_of_Gen_A;

static void *xmalloc (int size, const char *string)
{
  void *pointer;

  if( (pointer = malloc(size)) == NULL)
    {
      perror (string);
      exit (2);
    }

  return pointer;
}

static int krit_kon_mat (matrix_TYP **A, matrix_TYP **B, int anz,
			 matrix_TYP ***M, int *dim)
{
  matrix_TYP **AA, **BA, **AB,**Prod, *mat;
  int ab, aa, ba, pp,i,j;
  
  AB = solve_endo (A, B, anz, &ab);
  BA = solve_endo (B, A, anz, &ba);
  *dim = ab;
  *M = AB;
  if (ab == ba)
    { 
      AA = solve_endo (A, A, anz, &aa);
      if (ab == aa)
	{
	  Prod = (matrix_TYP **) xmalloc (aa * aa * sizeof (matrix_TYP *), "krit_kon_mat"); /* Prod = (matrix_TYP **) calloc (aa*aa, sizeof (matrix_TYP *)); */
	  for (i=0; i<aa; i++)
	    for (j=0; j<aa; j++)
	      Prod[i*aa + j] = mat_mul (AB[i], BA[j]);

	  mat = mat_to_line (Prod, aa*aa);

	  for (i=0; i<aa*aa; i++)
	    free_mat (Prod[i]);
	  free (Prod);

	  pp = tgauss (mat);		
	  free_mat (mat);

	  if (pp == aa)
	    {
	      if (ab != 0)
		{
		  for (i=0; i<ba; i++)
		    free_mat (BA[i]);
		  free (BA);
		}
	      if (aa != 0)
		{
		  for (i=0; i<aa; i++)
		    free_mat (AA[i]);
		  free (AA); 
		}
	      return 1;
	    }
	  else 
	    {
	      if (ab != 0)
		{
		  for (i=0; i<ba; i++)
		    free_mat (BA[i]);
		  free (BA);
		}
	      if (aa != 0)
		{
		  for(i=0; i<aa; i++)
		    free_mat (AA[i]);
		  free (AA); 
		} 
	      return 0;
	    }
	}
      else 
	{
	  if (ab != 0)
	    {
	      for (i=0; i<ba; i++)
		free_mat (BA[i]);
	      free (BA);
	    }
	  if (aa != 0)
	    { 
	      for (i=0; i<aa; i++)
		free_mat (AA[i]);
	      free (AA); 
	    }
	  return 0;
	}
    }
  else 
    {
      if (ab != 0)
	{
	  for (i=0; i<ba; i++)
	    free_mat (BA[i]);
	  free (BA);
	}
      if (aa != 0)
	{
	  for (i=0; i<aa; i++)
	    free_mat (AA[i]);
	  free (AA); 
	}
      return 0; 
    }
  
}


static matrix_TYP *suche_max_rank (matrix_TYP *mat_best, matrix_TYP **solve,
				   int dim, int step, int *found,
				   int rank_mat_best, int *rank)
{ 
  matrix_TYP **alternative, *m, *n, *B, *gau, *mat;
  int i, l, rank_m, rank_n, rank_B, rank_mat, z = 0;
  
  B = solve[step];
  gau = ggauss(B);
  rank_B = gau->rows;
  alternative = (matrix_TYP **) calloc (rank_B,sizeof(matrix_TYP *));
  free_mat (gau);

  for(l=1; l<rank_B+1  &&  rank_mat_best != B->rows; l++)
    {
      mat = imat_add (mat_best,B,1,l);		
      gau = ggauss (mat);
      rank_mat = gau->rows;
      free_mat (gau);
      if (rank_mat > rank_mat_best)
	{
	  rank_mat_best=rank_mat;
	  for (i=0; i<z; i++)
	    free_mat(alternative[i]);
	  z = 0;
	  if(mat_best != solve[0])
	    free_mat(mat_best);
	  mat_best=mat;
	}
      else
	{
	  if(rank_mat == rank_mat_best)
	    {
	      alternative[z] = mat;
	      z++;
	    }
	  else	
	    free_mat(mat);
	}	
    }
  *rank = rank_mat_best;
  if (rank_mat_best == B->rows)
    {
      *found = 1;
      for(i=0; i<z; i++)
	free_mat (alternative[i]);
      free (alternative);

      if (mat_best != solve[0])
	return mat_best;
      else
	return copy_mat(mat_best);
    }
  if(step == dim-1)
    {
      *found = 0;
      for (i=0; i<z; i++)
	free_mat (alternative[i]);
      free (alternative);
      return mat_best;
    }
  m = suche_max_rank (mat_best, solve, dim, step+1, found,
		      rank_mat_best, &rank_m);
  while (z > 0   &&   *found == 0)
    {
      n = suche_max_rank (alternative[z-1], solve, dim, step+1, found,
			  rank_mat_best, &rank_n);
      z--;
      if(rank_n > rank_m)
	{
	  free_mat (m);
	  m = n;
	  rank_m = rank_n;
	}
      else
	free_mat (n);
    }
  *rank = rank_m;

  for (i=0; i<z; i++)
    free_mat (alternative[i]);
  free (alternative);
  
  return m;
}


static matrix_TYP *suche_kon_mat (matrix_TYP **solve, int dim)
{
  int rank_mat_best, found = 0;
  matrix_TYP *A, *gau, *mat_best;
  
  A = solve[0];
  gau = ggauss (A);
  /* rank_mat_best=gau->rows; */
  if(dim > 1)
    mat_best = suche_max_rank (A, solve, dim, 1, &found, gau->rows,
			       &rank_mat_best);
  else 
    {
      rank_mat_best = gau->rows;
      mat_best = copy_mat (A);
    }
  free_mat (gau);
  if(rank_mat_best != A->rows)
    printf("Keine invertierbare Loesung gefunden.\n");
  return(mat_best);
}

/* Hier wird die Ordnung von m ausgerechnet, und zurueckgegeben.   */

static int order (matrix_TYP *m)
{
  int i = 1;
  matrix_TYP *n, *E;
  
  n = copy_mat (m);
  E = einheitsmatrix(m->rows);

  while(cmp_mat(n, E) != 0)
    {
      n = mat_muleq(n, m);
      
      i++;
    }

  free_mat(E);
  
  free_mat(n);
  
  return i;
}

static int **calc_traces_of_prod (bravais_TYP *G)
{
  int **mat, i, j;
  
  matrix_TYP *waste;
  
  mat = (int **) xmalloc(G->gen_no*sizeof(int *), "calc_traces_of_prod");
  
  for (i=0; i<G->gen_no; i++)
    mat[i] = (int *) xmalloc(G->gen_no*sizeof(int),"calc_traces_of_prod");

  /* Keine Lust Symmetrie der Matrix zu benutzen! */
  for (i=0; i<G->gen_no; i++)
    for (j=0; j<G->gen_no; j++)
      {
	waste = mat_mul(G->gen[i],G->gen[j]);
	
	mat[i][j] = trace(waste);
	
	free_mat(waste);
      }
  
  
  return mat;
}

static int *calc_traces_of_inv (bravais_TYP *G)
{
  int *mat, i;
  
  matrix_TYP *waste;
  
  mat = (int *) xmalloc(G->gen_no*sizeof(int), "calc_traces_of_inv");
  

  for (i=0; i<G->gen_no; i++)
      {
	waste = mat_inv(G->gen[i]);
	
	mat[i] = trace(waste);
	
	free_mat(waste);
      }
  
  
  return mat;
}

/* Get Conjugation Classes of "G" under Group "H". Store the con_class
   of each element in "list[]" and an "representant_of_con_class" in
   the array. Return_value is the last array. */

static int *con_classes (matrix_TYP **G, int order_of_G,
			 bravais_TYP *H, int *number_of_con_classes)
{
  matrix_TYP **orbit;
  int i, j, pos_in_list, length, *representant_of_con_class, *list;
  int orbit_options[6] = {4, 0, 0, 0, 0, 0};
  /* acting.from = Left,
     acting.by = Conjugation;*/

  *number_of_con_classes = 0;  

  list = (int *)xmalloc(order_of_G * sizeof(int), "con_classes");
  
  representant_of_con_class = (int *)xmalloc(order_of_G * sizeof(int),
					     "con_classes");
  /* Worst-Case-Allocation: Every element its own con_class. */
  
  for (i=0; i<order_of_G; i++)
    list[i] = -1;
  
  for (i=0; i<order_of_G; i++)
    if(list[i] == -1)
      {
	orbit = orbit_alg(G[i], H, NULL, orbit_options,
			  &length);

	representant_of_con_class[*number_of_con_classes] = i;
	
	for (j=0; j<length; j++)
	  {
	    pos_in_list = mat_search (orbit[j], G, order_of_G, mat_comp);
	    list[pos_in_list] = *number_of_con_classes;
	    
	    free_mat (orbit[j]);
	  }
	
	free (orbit);
	
	(*number_of_con_classes) ++;
	
      }
  
  representant_of_con_class = (int *) realloc(representant_of_con_class, *number_of_con_classes * sizeof(int));
  /* Just for the fun of saving a bit of memory. */

  free(list);

  return representant_of_con_class;
}



/* compute the trace of every element in "group" and store it in "trace_list". */

static int *compute_trace (matrix_TYP **group, int group_order)
{
  int i;
  
  int *trace_list = (int *) xmalloc (group_order * sizeof(int), "compute_trace");
  
  for (i=0; i<group_order; i++)
    trace_list[i] = trace (group[i]);
	
  return trace_list;
  
} 



/* compute the order of every element in "group" and store it in "order_list". */

static int *compute_order  (matrix_TYP **group, int group_order)
{
  int i;
  
  int *order_list = (int *) xmalloc (group_order * sizeof(int), "compute_order");
  
  for (i=0; i<group_order; i++)
    order_list[i] = order (group[i]);

  return order_list;
  
} 


static int end_test (matrix_TYP **matrix, matrix_TYP **A, matrix_TYP **B,
		     int number_of_el_in_A)
{
  matrix_TYP **AA;/*, *mat;*/
  int found, dim, i;
  
  found = krit_kon_mat(A, B, number_of_el_in_A, &AA, &dim);
  
  if (found == TRUE) 
    {
      /*mat = suche_kon_mat(AA,dim);	
       *matrix = mat;*/
      *matrix = suche_kon_mat(AA, dim);
      
      for(i=0; i<dim; i++)
	free_mat (AA[i]);
      free (AA);
    }
  
  return found;
}

/* returns centralizer of the group generated by elements "H" in
   "list_of_H". "list_of_H" is "number_of_H"-elements long. */

static bravais_TYP *get_centralizer(matrix_TYP **H, bravais_TYP *G,
				    int list_of_H[], int number_of_H)
{
  bravais_TYP *C1, *C2;

  matrix_TYP **waste;
  
  int orbit_options[6] = {4, 0, 3, 0, 0, 0};
  /*acting.and_calc_stab = TRUE;*/
  /* act.from = Left,
     act.by = Conjugation;*/
  int i, j, laenge;

  C1 = copy_bravais (G);
  
  C2 = init_bravais (G->dim);
  
  /* Dies laeuft in etwa darauf hinaus die Stabilisatoren aller "H"s zu
     schneiden. */
  for (i=0; i<number_of_H; i++)
    {
      waste = orbit_alg(H[list_of_H[i]], C1, C2, orbit_options, &laenge);

      for(j=0; j<laenge; j++)
	free_mat(waste[j]);
      free(waste);
    
      free_bravais(C1);
      
      C1 = C2;
      
      C2 = init_bravais (G->dim);

    }

  free_bravais (C2);
  
  return C1;
}


static int traces_are_ok(matrix_TYP *new_element, bravais_TYP *Group_A, matrix_TYP **Group_B, int next_gen)
{
  matrix_TYP *waste;

  int i; 
 
  waste = mat_inv(new_element);
  
  if (trace(waste) != traces_of_inv_A[next_gen])
    {
      free_mat (waste);
      
      return FALSE;
    }

  free_mat (waste);

  waste = mat_mul(new_element, new_element);
  
  if (trace(waste) != traces_of_prod_A[next_gen][next_gen])
    {
      free_mat (waste);
      return FALSE;
    }
      
  free_mat (waste);

  for (i=0; i<next_gen; i++)
    {
      waste = mat_mul(new_element, Group_B[image_of_Gen_A[i]]);
      
      if (trace(waste) != traces_of_prod_A[next_gen][i])
	{
	  free_mat (waste);
	  return FALSE;
	}
      
      free_mat (waste);
    }
  

  return TRUE;
  
}


/* function tests if mapping of the generators "Gen_A" in "Group_B" exists
   which projects "Gen_A[0],..,Gen_A[next_gen -1]" on
   "Group_B[image_of_Gen_A[0]],...,Group_B[image_of_Gen_A]". writes a
   possible result back into  the rest of the global array "image_of_Gen_A[]".
   Returns "1" if it succeeds and "0" if it does not. */
static int construct_Image_Gen_A (int next_gen, bravais_TYP *Gen_A,
				  int group_order, bravais_TYP *Gen_B,
				  matrix_TYP **Group_B)
{
  int *rep_con_class, number_of_con_classes,
    trace_next_gen, order_next_gen;
  
  matrix_TYP *matrix, **Bild;

  int i, j;
  
  bravais_TYP *centralizer;

  trace_next_gen = trace(Gen_A->gen[next_gen]);
  order_next_gen = order(Gen_A->gen[next_gen]);
  
  /* Stabilizer of trivial group is whole group. Maybe faster. */
  if (next_gen> 0)
    centralizer = get_centralizer(Group_B, Gen_B, image_of_Gen_A, next_gen);
  else
    centralizer = copy_bravais(Gen_B);
  
  rep_con_class = con_classes (Group_B, group_order, centralizer,
			       &number_of_con_classes);
  


  Bild = (matrix_TYP **)xmalloc( (next_gen+1) * sizeof(matrix_TYP *),"construct_Image_Gen_A");
  
  
  for (i=0; i<next_gen; i++)
    Bild[i] = Group_B[image_of_Gen_A[i]];

  
  for (i=0; i<number_of_con_classes; i++)
    if(trace_list_B[rep_con_class[i]] == trace_next_gen  &&
       order_list_B[rep_con_class[i]] == order_next_gen  &&
       traces_are_ok(Group_B[rep_con_class[i]], Gen_A, Group_B, next_gen) == TRUE)
      {
	Bild[next_gen] = Group_B[rep_con_class[i]];

	if(end_test (&matrix, Bild, Gen_A->gen, next_gen + 1) == TRUE)
	  {
	    free_mat(matrix);
	    
	    image_of_Gen_A[next_gen] = rep_con_class[i];
	    
	    if(next_gen + 1 == Gen_A->gen_no)
	      {
		free (Bild);
		free (rep_con_class);
		
		free_bravais (centralizer);
		
		return TRUE;
	      }
	    else
	      if (construct_Image_Gen_A (next_gen + 1, Gen_A, group_order, Gen_B, Group_B) == TRUE)
		{
		  free (Bild);
		  free (rep_con_class);
		  
		  free_bravais (centralizer);
		  
		  return TRUE;
		}

	  }

      }

  free (Bild);
  free (rep_con_class);
  
  free_bravais (centralizer);

  return FALSE;
}


matrix_TYP *suche_kand (bravais_TYP *Gen_A, bravais_TYP *Gen_B)
{
  matrix_TYP **Group_A, **Group_B;
  int group_order,laenge, i;
  
  matrix_TYP *matrix = NULL, **Bild;

  int orbit_options[6] = {0, 0, 0, 0, 0, 0};
  
  /*  acting.from = Left,
      acting.by = Normal_Operation,
      acting.on = Matrix_M,
      acting.until_orbit_length = 0,
      acting.and_calc_stab = FALSE,
      acting.until_stab_length = 0;
      acting.with_inverse = FALSE;
  */ 


/*  Check, if the groups have same dimension. */
 if (Gen_A->dim != Gen_B->dim)
    {
      fprintf (stdout, "Groups have different dimension.\n");
      
      return NULL;
    }

  /* Compute all elements of the two groups */
  Group_A = orbit_alg (Gen_A->gen[0], Gen_A, NULL, orbit_options, &group_order);

  Group_B = orbit_alg (Gen_B->gen[0], Gen_B, NULL, orbit_options, &laenge);

 /* We were only interested in the order of Group A,
    now the Generators are enough. */
  for (i=0; i<group_order; i++)
    free_mat (Group_A[i]);
  free(Group_A);

  /*  Check, if the groups have same order. */
  if (group_order != laenge)
    {
      fprintf (stdout, "Groups are of different order.\n");
      
      for (i=0; i<laenge; i++)
	free_mat (Group_B[i]);
      
      free (Group_B);
      
      return NULL;
    }

  /* This is necessary because of a "search_mat", at least this is
     what I believe. */
  mat_quicksort (Group_B, 0, group_order - 1, mat_comp);
  

  
  /*  orbit_options[0] = 4;*/
  /* acting.from = Left,
     acting.by = Conjugation;*/
  

  
  trace_list_B = compute_trace (Group_B, group_order);
  
  order_list_B = compute_order (Group_B, group_order);
  
 
  traces_of_prod_A = calc_traces_of_prod (Gen_A);
  traces_of_inv_A = calc_traces_of_inv (Gen_A);

  
  image_of_Gen_A = (int *) xmalloc (Gen_A->gen_no * sizeof(int), "suche_kand");


  if (construct_Image_Gen_A( 0, Gen_A,group_order, Gen_B, Group_B) == TRUE)
    {
      Bild = (matrix_TYP **)xmalloc(Gen_A->gen_no * sizeof(matrix_TYP *),"construct_Image_Gen_A");
      
      for (i=0; i<Gen_A->gen_no; i++)
	Bild[i] = Group_B[image_of_Gen_A[i]];
      
      
      if(end_test (&matrix, Bild, Gen_A->gen, Gen_A->gen_no) != TRUE)
	matrix = NULL;
      
      free(Bild);
    }
 
  
  free (image_of_Gen_A);
  
  for (i=0; i<group_order; i++)
    free_mat (Group_B[i]);
  free (Group_B);
  
  free(trace_list_B);
  
  free(order_list_B);
  
  
  for (i=0; i<Gen_A->gen_no; i++)
    free (traces_of_prod_A[i]);
  free (traces_of_prod_A);
  
 free (traces_of_inv_A);


 return matrix;
}




