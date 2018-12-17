
#include <signal.h>
#include "typedef.h"
#include "getput.h"
#include "tools.h"
#include "bravais.h"
#include "orbit.h"
#include "base.h"
#include "matrix.h"
#include "longtools.h"
#include "voronoi.h"
#include "sort.h"
#include "reduction.h"
#include "autgrp.h"
#include "ZZ_P.h"
#include "ZZ_irr_const_P.h"
#include "ZZ_lll_P.h"
#include "ZZ_cen_fun_P.h"
#include "ZZ_zclass_P.h"
#include "datei.h"

extern int INFO_LEVEL;


static int position(matrix_TYP **a,matrix_TYP *x,int n)
/* returns the first index i<n such that a[i] == x, and
   -1 if none exists */
{
  int i=0;

  while (i<n){
    if (mat_comp(a[i],x) == 0){
       return i;
    }
    i++;
  }
  return -1;
}


/*------------------------------------------------------------------------------- */
static matrix_TYP *is_conjugated_ZZ(n,new)
     ZZ_node_t *n, *new;
{
	int i,
	    Nanz;

	matrix_TYP *Tmp1,
		   *Tmp2,
		   *A,
		   *X = NULL,
	          **BASE,
                  **N;

	bravais_TYP *H=NULL;

	/* save the normailzer of n->colgroup in case it exits */
	N = n->col_group->normal;
	Nanz = n->col_group->normal_no;
	n->col_group->normal = NULL;
	n->col_group->normal_no = 0;

	/* bare in mind that is_z_equivalent has got an sideefect in
	calculating the normalizer of the bravais group of the first
	argument, and that it only calculates an isometry of formspaces. */
	Tmp1 = is_z_equivalent_datei(n->col_group,
	n->group,new->col_group,new->group,&n->perfect,&n->perfect_no);

	/* look whether we already got the bravais group n->brav */
	if (n->brav == NULL){
		n->brav = bravais_group(n->col_group,TRUE);
	}

	if (n->col_group->normal != NULL &&
	    n->col_group->normal_no > 0){
		/* stick the normalizer to the bravais group */
		n->brav->normal = n->col_group->normal;
		n->brav->normal_no = n->col_group->normal_no;
	}

	n->col_group->normal = N;
	n->col_group->normal_no = Nanz;

	if (Tmp1 != NULL){
		/* we should get a stabchain for n->col_group */
		if (n->stab_chain == NULL){
			BASE = get_base(n->col_group);
			n->stab_chain=strong_generators(BASE,n->col_group,FALSE);
			/* free the base again */
			for (i=0;i<n->col_group->dim;i++){
				free_mat(BASE[i]);
			}
			free(BASE);
		}

		Tmp2 = long_mat_inv(Tmp1);
		H = konj_bravais(new->col_group,Tmp2);
		free_mat(Tmp2);

		/* now we are in the position to look if H and n->col_group
		are conjugated */
		Nanz = n->brav->gen_no + n->brav->normal_no;
		N = (matrix_TYP **) malloc(Nanz * sizeof(matrix_TYP *));
		for (i=0;i<Nanz;i++){
			if (i<n->brav->gen_no){
				N[i] = n->brav->gen[i];
			}
			else{
				N[i] = n->brav->normal[i-n->brav->gen_no];
			}
		}

		if (n->brav->normal_no == 0){
			fprintf(stderr,
			        "ZZ_zclass: bravais group has no normalizer\n");
			exit(3);
		}

		A = conjugated(n->col_group,H,N,Nanz,
				n->stab_chain);

		free(N);
		
		if (A == NULL){
		   free_mat(Tmp1);
		}
		else{
		   X = mat_mul(Tmp1, A);
		   free_mat(A);
		   free_mat(Tmp1);
		}
	}

	if (H!=NULL){
		free_bravais(H);
	}


	return (X);
}



/*------------------------------------------------------------------------------- */
int deal_with_ZCLASS(data, tree, father, new)
     ZZ_data_t *data;
     ZZ_tree_t *tree;
     ZZ_node_t *father, *new;
{

	int f,
	    g;

	bravais_TYP *H;

	ZZ_node_t *n;
	
	matrix_TYP *X;
	


	/* now calculate the (integral) representation
	on the new lattice (bare in mind that it is
	row invariant. There is a flaw: normalizer and
	centralizer won't be correct */
	f = tree->root->group->normal_no;
	tree->root->group->normal_no=0;
	g = tree->root->group->cen_no;
	tree->root->group->cen_no=0;
	new->group = konj_bravais(tree->root->group,new->U);
	tree->root->group->normal_no = f;
	tree->root->group->cen_no = g;

	/* the second flaw of konj_bravais is the formspace */
	for (f=0;f<new->group->form_no;f++)
		new->group->form[f]->kgv = 1;
	long_rein_formspace(new->group->form,new->group->form_no,1);

	new->col_group = tr_bravais(new->group,1,FALSE);

	/* see if we've got an Z-equivalent rep. already */
	n = tree->root;
	while (n != NULL){
	        X = is_conjugated_ZZ(n,new);
		if (X != NULL){
		        free_mat(X);
			return TRUE;
		}
		n = n->next;
	}
			
	/* exit(3); */

	return FALSE;
}




/* -------------------------------------------------------------------------- */
matrix_TYP *special_deal_with_zclass(ZZ_tree_t *tree,
                                     bravais_TYP *group,
                                     int *nr)
{
   ZZ_node_t *n,
             *new;

   matrix_TYP *X;


   nr[0] = 0;
   new = (ZZ_node_t *)calloc(1, sizeof(ZZ_node_t));
   new->col_group = group;
   new->group = tr_bravais(group, 1, FALSE);

   n = tree->root;
   while (n != NULL){
      X = is_conjugated_ZZ(n, new);
      if (X != NULL){
         free_bravais(new->group);
         new->col_group = NULL;
         free(new);
         n = NULL;
         return(X);
      }
      nr[0]++;
      n = n->next;
   }

   fprintf(stderr,"ERROR in special_deal_with_zclass\n");
   exit(4);
}




/*------------------------------------------------------------------------------- */
static int is_contained(matrix_TYP *U,matrix_TYP *V)
{

	int i,
	    j,
	    k,
	    sum;

	for (i=0;i<U->cols;i++){
		for (j=0;j<U->rows;j++){
			sum = 0;
			for (k=0;k<U->rows;k++){
				sum += (V->array.SZ[i][k] * U->array.SZ[k][j]);
			}
			if ((sum % V->kgv) != 0)
				return FALSE;
		}
	}

	return TRUE;
/*
	matrix_TYP *tmp;

	tmp = mat_mul(V,U);

	Check_mat(tmp);

	if (tmp->kgv == 1){
		free_mat(tmp);
		return TRUE;
	}
	else{
		free_mat(tmp);
		return FALSE;
	}
*/

}



/*------------------------------------------------------------------------------- */
static int suche_mat(matrix_TYP *mat,
                     matrix_TYP **liste,
                     int anz)
{
   int i;

   for (i = 0; i < anz; i++){
      if (cmp_mat(mat, liste[i]) == 0)
         return(i);
   }
   return(-1);
}



/*------------------------------------------------------------------------------- */
int in_bahn(matrix_TYP *lattice,
            ZZ_node_t *father,
            int *i)
{
   int flag = 0;



   for (i[0] = 0; i[0] < father->N_no_orbits && !flag; i[0]++){
      flag = (mat_search(lattice, father->N_orbits[i[0]], father->N_lengths[i[0]], mat_comp) != -1);
   }

   return(flag);
}



/*------------------------------------------------------------------------------- */
int orbit_under_normalizer(data,tree,father,new,ii,jj,inzidenz,nr,nnn)
     ZZ_data_t *data;
     ZZ_tree_t *tree;
     ZZ_node_t *father, *new;
     int ii, jj;
     QtoZ_TYP *inzidenz;
     int *nr;
     ZZ_node_t **nnn;
{

	int i__,
	    j__,
	    i, pos, oflag = 0,
	    father_flag = FALSE,
	    flag = FALSE;

	static int orb[6];

	matrix_TYP *tmp;

	ZZ_node_t *p;
	
	ZZ_couple_t *laeufer;

	bravais_TYP *G;

	orb[0]=0;
	orb[1]=4;

	tmp = copy_mat(new->U);
	long_row_hnf(tmp);
	ZZ_transpose_array(tmp->array.SZ,tmp->cols);

	/* search the whole tree */
	p = tree->root;
	while (p!= NULL){
	   for (i__=0;i__<p->N_no_orbits && p->level==father->level && !flag;i__++){
	      flag = (mat_search(tmp,p->N_orbits[i__], p->N_lengths[i__],mat_comp) != -1);
	   }
	   if (flag) 
	      break;
	   p = p->next;
	}

	/* search the father
	p = father;
	for (i__=0;i__<p->N_no_orbits && !father_flag;i__++){
		father_flag = (mat_search(tmp,p->N_orbits[i__],
				p->N_lengths[i__]) != -1);
	} */
	free_mat(tmp);
	
	/* next 14 lines: oliver 10.8.00: graph for QtoZ */
	if (flag && GRAPH){
	   laeufer = p->child;
	   for (i = 0; i < p->N_no_orbits - i__; i++){
	      laeufer = laeufer->elder;
	      if (laeufer == NULL){
	         fprintf(stderr,"ERROR 1 in orbit_under_normalizer!\n");
	         exit(2);
	      }
	   }
	   nnn[0] = laeufer->he;
	   nr[0] = suche_mat(laeufer->he->U, inzidenz->gitter, inzidenz->anz);
	   if (nr[0] == -1){
	      fprintf(stderr,"ERROR 2 in orbit_under_normalizer!\n");
	      exit(3);
	   }
	}
	/* else { */
	p = father;
	if (!flag){
	        /* p = father; */
		tmp = tr_pose(new->U);
		if (p->N_no_orbits == 0){
			p->N_orbits = (matrix_TYP ***) malloc(1 *
					sizeof(matrix_TYP **));
			p->N_lengths = (int *) malloc(1 * sizeof(int));
		}
		else{
			p->N_orbits = (matrix_TYP ***) realloc(p->N_orbits,
				(p->N_no_orbits+1)* sizeof(matrix_TYP **));
			p->N_lengths = (int *) realloc(p->N_lengths,
				(p->N_no_orbits + 1) * sizeof(int));
		}

		G = init_bravais(tree->root->group->dim);

		/* tilman: changed on 06.05.1998:
                G->gen_no = tree->root->col_group->normal_no;
		G->gen = tree->root->col_group->normal; */		
		G->gen_no = tree->root->group->normal_no;
		G->gen = tree->root->group->normal;

		p->N_orbits[p->N_no_orbits] = orbit_alg(tmp,G,NULL,orb,
				&p->N_lengths[p->N_no_orbits]);

		/* sort it */
		mat_quicksort(p->N_orbits[p->N_no_orbits],0,
				p->N_lengths[p->N_no_orbits]-1,mat_comp);

		p->N_no_orbits++;

		G->gen_no = 0;
		G->gen = NULL;
		free_bravais(G);
		free_mat(tmp);
	}

	if (!flag){
		fprintf(stderr,"father: level %d number of orbits %d lengths %d\n",
			p->level,p->N_no_orbits,p->N_lengths[p->N_no_orbits-1]);
	}

	if (INFO_LEVEL & 8){
		if (flag || father_flag){
			fprintf(stderr,"returned TRUE\n");
		}
		else{
			fprintf(stderr,"returned FALSE\n");
		}
	}
	
	return (flag || father_flag);

}



/* ---------------------------------------------------------------------------------------- */
matrix_TYP *konjugierende(matrix_TYP *Li,
                          bravais_TYP *G,
                          ZZ_node_t *n)
{
   ZZ_node_t *new;

   matrix_TYP *X;

   bravais_TYP *group;

   int f, g;



   f = G->normal_no; G->normal_no = 0;
   g = G->cen_no; G->cen_no = 0;
   group = konj_bravais(G, Li);
   G->normal_no = f;
   G->cen_no = g;
   for (f = 0; f <  group->form_no; f++)
      group->form[f]->kgv = 1;
   long_rein_formspace(group->form, group->form_no, 1);

   new = (ZZ_node_t *)calloc(1, sizeof(ZZ_node_t));
   new->col_group = group;
   new->group = tr_bravais(group, 1, FALSE);

   X = is_conjugated_ZZ(n, new);

   /* clean */
   free_bravais(new->group);
   free_bravais(new->col_group);
   free(new);

   return(X);
}











