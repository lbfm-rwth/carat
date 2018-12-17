#include "typedef.h"
#include "voronoi.h"
#include "reduction.h"
#include "autgrp.h"
#include "bravais.h"
#include "symm.h"
#include "getput.h"
#include "matrix.h"
#include "longtools.h"
#include "polyeder.h"
#include "tools.h"
#include "sort.h"
#include "longtools.h"
#include "datei.h"

extern int INFO_LEVEL;

int max_diagonal_entry(matrix_TYP *A)
/* returns the entry |A[i][i]| such that |A[i][i]| >= |A[j][j]| for all j */
{
   int max=0,
       i;

   if (A->cols != A->rows){
      printf("error in max_diagonal_entry\n");
      exit(3);
   }

   for (i=0;i<A->cols;i++){
      if ((max==0) || (max < abs(A->array.SZ[i][i]))){
         max = abs(A->array.SZ[i][i]);
      }
   }

   return max;
}

matrix_TYP *extends_to_isometry(
matrix_TYP **hforms,matrix_TYP *HSV,int anz_hneighbours,
matrix_TYP **gforms,matrix_TYP *GSV,int anz_gneighbours,
int fdim,int offset)
/* the situation is as follows:
   the matrices hform and gform are isometric, now looks whther
   we can find an isometry which transforms the set hneighbours
   do the set gneighbours */
{
   matrix_TYP *erg=NULL,
              *tmp;

   int i=offset;

   if (INFO_LEVEL & 12){
      fprintf(stderr,"entering extends_to_isometry\n");
   }
   if (anz_hneighbours != anz_gneighbours){
      fprintf(stderr,"WARNING IN extends_to_isometry: different number\n");
      fprintf(stderr,"                                of forms.\n");
   }

   while ((erg==NULL) && (i<anz_hneighbours)){
     /* swap the i-th entry of gneighbours to the offset-th entry */
     tmp = gforms[offset];
     gforms[offset] = gforms[i];
     gforms[i] = tmp;

     /* test whether this setting is possible at all */
     tmp = isometry(hforms,gforms,offset+1,HSV,GSV,NULL,0,NULL);

     if (tmp!=NULL){
       /* changed 16.10.01 tilman: the first condition is wrong, because
	  it uses properties of the direction (i.e. to be a base of the
	  space of invariant form) which they do not fullfill */
       /* if ((offset<(fdim-1)) && (offset < (anz_hneighbours - 2))){*/
       if (offset < (anz_hneighbours - 2)){
           /* do the recursion */
           free_mat(tmp);
           erg = extends_to_isometry(hforms,HSV,anz_hneighbours,
                                  gforms,GSV,anz_gneighbours,fdim,offset+1);
        }
        else{
           erg = tmp;
        }
     }

     /* swap back */
     tmp = gforms[i];
     gforms[i] = gforms[offset];
     gforms[offset] = tmp;
     i++;
   }

   return erg;
}


int neighbours(matrix_TYP ***perf,bravais_TYP *G,matrix_TYP **Ftr,
               matrix_TYP *tr_bifo,matrix_TYP *SV,int min)
{
  int i,j,k,l, Fdim, Gdim;
  int lc, rc, g,
      min_x;
  int Pmin, Vanz, Wanz;

  matrix_TYP **V,
              *N,
              *P,      /* the G-perfect form */
              *X;
  polyeder_TYP *pol;
  wall_TYP **W;
  int *ww, *pw;

  /* initial settings */
  P = perf[0][0];

  Fdim = G->form_no;
  Gdim = G->form[0]->cols;
  /**********************************************************************\
  | calculate the vertices of the Voronoidomain in F(G^{tr})
  | and the forms in F(G) corresponding to the walls of the Voronoidomain.
  | The coordinates of these forms with respect to the basis given in G->form
  | are obtained as the coordinates of the vetrices of the polyeder_TYP *pol.
  \**********************************************************************/
  V = voronoi_vertices(P, G, &Vanz, &Pmin, &i);
  Wanz = Vanz;
  if( (W = (wall_TYP **)malloc(Wanz *sizeof(wall_TYP *))) == NULL)
  {
     printf("malloc of W failed in neighbours\n");
     exit(2);
  }
  if( (ww = (int *)malloc(Fdim *sizeof(int))) == NULL)
  {
     printf("malloc of ww failed in neighbours\n");
     exit(2);
  }
  if( (pw = (int *)malloc(Fdim *sizeof(int))) == NULL)
  {
     printf("malloc of pw failed in neighbours\n");
     exit(2);
  }
  for(i=0;i<Vanz;i++)
  {
    W[i] = init_wall(Fdim);
    form_to_vec(ww, V[i], Ftr, Fdim, &k);
    for(j=0;j<Fdim; j++)
      for(k=0;k<Fdim; k++)
         W[i]->gl[j] += ww[k] * tr_bifo->array.SZ[j][k];
    normal_wall(W[i]);
    free_mat(V[i]);
  }
  free(V);
  pol = first_polyeder(W, Wanz);
  if(pol == NULL)
  {
     printf("Error in neighbours:  P not G-perfekt\n");
     exit(3);
  }
  for(i=0;i<Wanz;i++)
  {
    refine_polyeder(pol, W[i]);
    free_wall(&W[i]);
  }
  free(W);

  /*
  put_polyeder(pol);
  */

  /******************************************************************\
  | Each G-perfect form A, that is a neighbour of P, is (up to
  | multiplication with a positive scalar) given by
  | N = lc *P + rc * X,
  | where X is a form defined by a vertex of the polyeder pol.
  | The coefficients are calculated with the function 'voronoi_neighbour'. 
  \******************************************************************/
  form_to_vec(pw, P, G->form, Fdim, &k);

  /* getting memory for the result */
  perf[0] = (matrix_TYP **) realloc(perf[0],(pol->vert_no + 1)
                                            * sizeof(matrix_TYP *));

  X = init_mat(Gdim, Gdim, "");
  X->flags.Integral = X->flags.Symmetric = TRUE;
  X->flags.Diagonal = FALSE;

  for(i=0;i<pol->vert_no;i++)
  {
    /******************************************************\
    | Calculate X
    \******************************************************/
    for(j=0;j<Gdim;j++)
      for(k=0;k<=j;k++)
      {
        X->array.SZ[j][k] = 0;
        for(l=0;l<Fdim;l++)
          X->array.SZ[j][k] +=pol->vert[i]->v[l] * G->form[l]->array.SZ[j][k];
        X->array.SZ[k][j] = X->array.SZ[j][k];
      }
    N = voronoi_neighbour(P, X, Pmin, &lc, &rc);

    if(N != NULL){
       perf[0][i+1] = N;
    }
    else{
       /* we got a blind direction here */
       /* so the neighbour in this sense is the intersection of the
          boundary of the cone with the line through P with slope
          X */
       N = scal_pr(SV,X,TRUE);
       min_x = N->array.SZ[0][0];
       for (j=0;j<N->rows;j++){
          if (N->array.SZ[j][j] < min_x) min_x = N->array.SZ[j][j];
       }
       perf[0][i+1] = imat_add(P,X,min_x,min);
       free_mat(N);
    }

    /* divide this result by the gcd of all entries */
    g = perf[0][i+1]->array.SZ[0][0];
    for (j=0;j<perf[0][i+1]->rows;j++){
       for (k=0;k<perf[0][i+1]->cols;k++){
          g = GGT(perf[0][i+1]->array.SZ[j][k],g);
       }
    }
    for (j=0;j<perf[0][i+1]->rows;j++){
       for (k=0;k<perf[0][i+1]->cols;k++){
          perf[0][i+1]->array.SZ[j][k] = perf[0][i+1]->array.SZ[j][k]/g;
       }
    }
  }

  free_mat(X);
  free(pw);
  free(ww);
/*********************************
  put_polyeder(pol);
***********************************/

  /* changed 14.04.99 to save this variable from being freed, it
     causes trouble on some machines */
  i = pol->vert_no;
  free_polyeder(pol);

  /* return the number of neighbours */
  return i;

}

void transform_pair(bravais_TYP *H,bravais_TYP *Htr,matrix_TYP *x)
{
  int i;

  matrix_TYP *xi,
             *xitr,
             *xtr,
             *tmp;

  Check_mat(x);
  xi=long_mat_inv(x);
  xitr=tr_pose(xi);
  xtr=tr_pose(x);

  if (INFO_LEVEL & 4){
     put_mat(x,NULL,"x",2);
     put_mat(xtr,NULL,"xtr",2);
     put_mat(xi,NULL,"xi",2);
     put_mat(xitr,NULL,"xitr",2);
  }

  /* firstly transform H->gen */
  for (i=0;i<H->gen_no;i++){
     tmp = mat_mul(xitr,H->gen[i]);
     free_mat(H->gen[i]);
     H->gen[i] = mat_mul(tmp,xtr);
     free_mat(tmp);
  }

  /* now transform Htr->gen */
  for (i=0;i<Htr->gen_no;i++){
     tmp = mat_mul(x,Htr->gen[i]);
     free_mat(Htr->gen[i]);
     Htr->gen[i] = mat_mul(tmp,xi);
     free_mat(tmp);
  }

  /* now transform H->form */
  for (i=0;i<H->form_no;i++){
     tmp = mat_mul(x,H->form[i]);
     free_mat(H->form[i]);
     H->form[i] = mat_mul(tmp,xtr);
     free_mat(tmp);
  }

  /* now transform Htr->form */
  for (i=0;i<Htr->form_no;i++){
     tmp = mat_mul(xitr,Htr->form[i]);
     free_mat(Htr->form[i]);
     Htr->form[i] = mat_mul(tmp,xi);
     free_mat(tmp);
  }

  /* old version
  tmp = konj_bravais(Htr,x);
  tmp = konj_bravais(H,xitr); */

  /* inserted: tilman 26.08.98 to get a better
     basis for the formspace of H,Htr */
  long_rein_formspace(H->form,H->form_no,1);
  long_rein_formspace(Htr->form,Htr->form_no,1);

  free_mat(xtr);
  free_mat(xi);
  free_mat(xitr);

  return;

}

matrix_TYP *is_z_equivalent(bravais_TYP *G,bravais_TYP *G_tr,
                            bravais_TYP *H,bravais_TYP *H_tr)
/* IMPORTANT: G,H HAVE TO BE BRAVAIS-GROUPS!!!!,
   AND G_tr, H_tr HAVE TO BE THE TRANSPOSED GROUPS OF G,H RESP.
   This function returns an integral matrix A such that G^A = H, if
   such a matrix exists. Otherwise NULL is returned */
{
   int i,
       j,
       dim,
       epsilon = 101,
       ming,
       minh,
       prime = 1949,
       anz_gperfect,
       anz_hneighbours,
       anz_gneighbours,
       max;          /* will hold the maximal diagonal entry of a form */

   voronoi_TYP **gp;       /* normalizer returns a list of all
                              perfect forms via a voronoi_TYP */

	bravais_TYP *GB;          /* the bravais group of G, for temporary purposes */

   matrix_TYP  *gtrbifo,     /* the tracebifo of G,G_tr */
               *htrbifo,     /* the tracebifo of H,H_tr */
               *tmp,
               *tmp2,
               *id,
               *gperfect,    /* will hold ONE G-perfect form modulo the
                              operation of the normalizer */
               *hperfect,    /* will hold ONE H-perfect form */
               *HSV,         /* will hold the short vectors of hperfect */
               *GSV,
              **hneighbours, /* holds all neighbours of hperfect */
              **gneighbours,
               *h_pair=NULL,      /* h will be transformed in a way such that
                                we get a "nice" basis of the formspace */
               *conj;        /* the matrix that conjugates these both groups,
                              if this exists */

   if (INFO_LEVEL & 4){
      fprintf(stderr,"Entered is_z_equivalent\n");
   }

   dim = G->dim;

   /* checking all the trivialities */
   if (G->dim != H->dim){
      conj = NULL;
      if (INFO_LEVEL & 4){
         fprintf(stderr,"the groups don't even live in the same universe\n");
      }
   }
   else if(G->form_no != H->form_no){
      conj = NULL;
      if (INFO_LEVEL & 4){
         fprintf(stderr,"the groups are not conjugated\n");
         fprintf(stderr,"Dimension of the Formspace of G: %d\n",G->form_no);
         fprintf(stderr,"Dimension of the Formspace of H: %d\n",H->form_no);
      }
   }
   /* the order is an criterion if it's known */
   else if((G->order != H->order) && (G->order !=0) && (H->order !=0)){
      conj = NULL;
      if (INFO_LEVEL & 4){
         fprintf(stderr,"the groups are not conjugated\n");
         fprintf(stderr,"Order of G: %d\n",G->order);
         fprintf(stderr,"Order of H: %d\n",H->order);
      }
   }
   else{
      /* now start the real bussiness */

      /* firstly calculate the trace_bifo for G,H */
      gtrbifo = trace_bifo(G->form,G_tr->form,G->form_no);
      htrbifo = trace_bifo(H->form,H_tr->form,H->form_no);

      /* inserted tilman 30.5.97 */
      for (i=0;i<G->form_no;i++){
         Check_mat(G->form[i]);
         Check_mat(H->form[i]);
      }

      /* output for debugging
      put_mat(gtrbifo,NULL,"gtrbifo",2);
      put_mat(htrbifo,NULL,"htrbifo",2); */

      /* the two trace bifos should have the same elementary devisors */
      tmp = long_elt_mat(NULL,gtrbifo,NULL);
      tmp2 = long_elt_mat(NULL,htrbifo,NULL);

      if (mat_comp(tmp,tmp2) == 0){
         free_mat(tmp);
         free_mat(tmp2);

         id = init_mat(dim,dim,"1");

         if (INFO_LEVEL & 4){
            fprintf(stderr,"\n");
         }

         /* now calculate one G-perfect form */
         tmp = rform(G->gen,G->gen_no,id,epsilon);
         gperfect = first_perfect(tmp,G,G_tr->form,gtrbifo,&ming);
         free_mat(tmp);

         /* now calculate one H-perfect form, and the short vectors
            of it */
         tmp = rform(H->gen,H->gen_no,id,epsilon);
         hperfect = first_perfect(tmp,H,H_tr->form,htrbifo,&minh);
         free_mat(tmp);

         /* reduce H by pair_reduction, and transform H_tr accordingly */
         h_pair = init_mat(H->dim,H->dim,"1");
         tmp = pair_red(hperfect,h_pair);
         free_mat(hperfect);

         transform_pair(H,H_tr,h_pair);

         free_mat(htrbifo);
         htrbifo = trace_bifo(H->form,H_tr->form,H->form_no);
         hperfect = first_perfect(tmp,H,H_tr->form,htrbifo,&minh);
         free_mat(tmp);

         max = max_diagonal_entry(hperfect);
         HSV = short_vectors(hperfect,max,0,0,0,&i);

         if (INFO_LEVEL & 4){
            fprintf(stderr,"Got perfect forms for G and H.\n");
         }

         /* output for debugging */
         if (INFO_LEVEL & 4){
            put_mat(hperfect,NULL,"hperfect",2);
            put_mat(htrbifo,NULL,"htrbifo",2);
            printf("maximaler diagonaleintrag von hperfect %d\n",max);
         }

         /* now calculate all perfect neigbours of hperfect */
         hneighbours = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
         hneighbours[0] = hperfect;
         anz_hneighbours = neighbours(&hneighbours,H,H_tr->form,htrbifo,
                                                                HSV,minh);

         if (INFO_LEVEL & 4){
            fprintf(stderr,"Calculated the neighbours of hperfect.\n");
         }

			GB = bravais_group(G,TRUE);
         /* and calculate all G-perfect forms which represent the orbit
            of the normalizer on them, i.e. the quotient graph */
         gp = normalizer(gperfect,GB,G_tr,prime,&anz_gperfect);
			G->normal = GB->normal; G->normal_no = GB->normal_no;
			GB->normal = NULL; GB->normal_no = 0;
			free_bravais(GB);

         if (INFO_LEVEL & 4){
            fprintf(stderr,"Calculated the normalizer of G.\n");
         }

         /* now search for a G-perfect form which is isometric to the
            H-perfect form and all of its neighbours are also */
         conj = NULL;
         i = 0;
         while ((i<anz_gperfect) && (conj == NULL)){

            /* calculate the short vectors of the G-perfect form in question */
            GSV = short_vectors(gp[i]->gram,max,0,0,0,&j);

            /* let's see whether there is an isometry between this form
               and hperfect */
            conj = isometry(&hperfect,&gp[i]->gram,1,HSV,GSV,NULL,0,NULL);

            if (INFO_LEVEL & 4){
               fprintf(stderr,"calculating isometries\n");
               fprintf(stderr,"%d-th of %d\n",i+1,anz_gperfect);
               if (conj!=NULL){
                  fprintf(stderr,"found a new isometry.\n");
                  put_mat(conj,NULL,"conj",2);
                  put_mat(hperfect,NULL,"hperfect",2);
                  put_mat(gp[i]->gram,NULL,"gp[i]->gram",2);
               }
            }

            /* output for debugging
            put_mat(hperfect,NULL,"hperfect",2);
            put_mat(gp[i]->gram,NULL,"gp[i]->gram",2);
            put_mat(conj,NULL,"conj",2);
            put_mat(HSV,NULL,"HSV",2);
            put_mat(GSV,NULL,"GSV",2);  */

            if (conj != NULL){
               /* there is an isometry, let's see whether it respects the
                  neighbours also */
               /* now calculate all perfect neigbours of hperfect */
               gneighbours = (matrix_TYP **) malloc(sizeof(matrix_TYP));
               gneighbours[0] = gp[i]->gram;
               anz_gneighbours = neighbours(&gneighbours,G,G_tr->form,gtrbifo,
                                                        GSV,ming);

               free_mat(conj);
               conj = NULL;

               if (anz_hneighbours == anz_gneighbours){

                  if (INFO_LEVEL & 4){
                     printf("anz_hneighbours %d\n",anz_hneighbours);
                     for (j=0;j<=anz_hneighbours;j++){
                        put_mat(hneighbours[j],NULL,"hneighbours[j]",2);
                        put_mat(gneighbours[j],NULL,"gneighbours[j]",2);
                     }
                  }

                  conj = extends_to_isometry(hneighbours,HSV,
                                         anz_hneighbours+1,gneighbours,GSV,
                                         anz_gneighbours+1,H->form_no,1);
  

               }
               /* cleaning up memory */
               for (j=0;j<anz_gneighbours;j++){
                  free_mat(gneighbours[j+1]);
               }
               free(gneighbours);

            }

            free_mat(GSV);
            i++;
         }

         /* clean up memory used only here */
         free_mat(id);
         free_mat(hperfect);
         free_mat(gperfect);
         free_mat(HSV);
         for (i=0;i<anz_gperfect;i++){
            clear_voronoi(gp[i]);
            free(gp[i]);
         }
         free(gp);
         for (i=0;i<anz_hneighbours;i++){
            free_mat(hneighbours[i+1]);
         }
         free(hneighbours);

      }
      else{
         /* the groups are not conjugated, so clean up the memory used
            so far and return NULL */
         if (INFO_LEVEL & 4){
            put_mat(tmp,NULL,"tmp",2);
            put_mat(tmp2,NULL,"tmp2",2);
            fprintf(stderr,
                    "The groups are not conjugated: elementary divisors\n");
         }
         free_mat(tmp);
         free_mat(tmp2);
         conj = NULL;
      }

      /* clean up memory used only in this stage */
      free_mat(gtrbifo);
      free_mat(htrbifo);
   }

   /* retransform H, H_tr to the original value, and fiddle around with
      conj a bit */
   if (h_pair != NULL){
      tmp = long_mat_inv(h_pair);
      transform_pair(H,H_tr,tmp);
      if (conj != NULL){
         mat_muleq(conj,h_pair);
      }
      free_mat(tmp);
      free_mat(h_pair);
   }

   /* we will return the transposed matrix if there was one at all */
   if (conj != NULL){
      tmp = tr_pose(conj);
      free_mat(conj);
      conj = tmp;
   }

   return conj;
}
