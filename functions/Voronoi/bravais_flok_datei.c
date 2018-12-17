#include "typedef.h"
#include "voronoi.h"
#include "getput.h"
#include "matrix.h"
#include "bravais.h"
#include "longtools.h"
#include "symm.h"
#include "autgrp.h"
#include "voronoi.h"
#include "reduction.h"
#include "sort.h"
#include "tools.h"
#include "datei.h"

extern int INFO_LEVEL;
 
matrix_TYP *is_z_equivalent_datei(bravais_TYP *G,bravais_TYP *G_tr,
          bravais_TYP *H,bravais_TYP *H_tr,voronoi_TYP ***gp,int *anz_gperfect)
{
   int i,
       j,
       dim,
       epsilon = 101,
       ming,
       minh,
       prime = 1949,
       anz_hneighbours,
       anz_gneighbours,
       max;          /* will hold the maximal diagonal entry of a form */

	bravais_TYP *GB;          /* the bravais group of G for temporary purposes */
	
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
      fprintf(stderr,"Entered is_z_equivalent_ZZ\n");
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
         if (gp[0] == NULL && anz_gperfect[0] == 0){
            tmp = rform(G->gen,G->gen_no,id,epsilon);
            gperfect = first_perfect(tmp,G,G_tr->form,gtrbifo,&ming);
            free_mat(tmp);
         }
         else{
            gperfect = NULL;
         }

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

         /* and calculate all G-perfect forms which represent the orbit
            of the normalizer on them, i.e. the quotient graph */
         if (gp[0]== NULL && anz_gperfect[0]==0){
				 GB = bravais_group(G,TRUE);
				 gp[0] = normalizer(gperfect,GB,G_tr,prime,anz_gperfect);
				 G->normal = GB->normal; G->normal_no = GB->normal_no;
				 GB->normal = NULL; GB->normal_no = 0;
				 free_bravais(GB);
			}

         if (INFO_LEVEL & 4){
            fprintf(stderr,"Calculated the normalizer of G.\n");
         }

         /* now search for a G-perfect form which is isometric to the
            H-perfect form and all of its neighbours are also */
         conj = NULL;
         i = 0;
         while ((i<anz_gperfect[0]) && (conj == NULL)){

            /* calculate the short vectors of the G-perfect form in question */
            GSV = short_vectors(gp[0][i]->gram,max,0,0,0,&j);

            /* let's see whether there is an isometry between this form
               and hperfect */
            conj = isometry(&hperfect,&gp[0][i]->gram,1,HSV,GSV,NULL,0,NULL);

            if (INFO_LEVEL & 4){
               fprintf(stderr,"calculating isometries\n");
               fprintf(stderr,"%d-th of %d\n",i+1,anz_gperfect[0]);
               if (conj!=NULL){
                  fprintf(stderr,"found a new isometry.\n");
                  put_mat(conj,NULL,"conj",2);
                  put_mat(hperfect,NULL,"hperfect",2);
                  put_mat(gp[0][i]->gram,NULL,"gp[i]->gram",2);
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
               /* because we have an isometry, the forms will have the same
                  minimum */
               ming = minh;

               /* now calculate all perfect neigbours of hperfect */
               gneighbours = (matrix_TYP **) malloc(sizeof(matrix_TYP));
               gneighbours[0] = gp[0][i]->gram;
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
         if (gperfect != NULL) free_mat(gperfect);
         free_mat(HSV);
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
