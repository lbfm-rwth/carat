/* author: Oliver Heidbuechel */
/* last change: 13.09.2000 */


#include <ZZ.h>
#include <typedef.h>
#include <bravais.h>
#include<getput.h>
#include<matrix.h>
#include<longtools.h>
#include<tools.h>
#include"zass.h"
#include"symm.h"
#include"base.h"
#include"autgrp.h"
#include"voronoi.h"
#include"polyeder.h"
#include <graph.h>


/* ---------------------------------------------------------------------------- */
/* Fill the structure Q_data_TYP for G, i.e. calculate various informations we  */
/* need for calculating the indlusiongraph of subgroup relations in a Q-class   */
/* pres has to be a correct presentation for G in Gl_n(Z)                       */
/* ---------------------------------------------------------------------------- */
Q_data_TYP *get_Q_data(bravais_TYP *G,
                       matrix_TYP *pres,
                       boolean l_option)
{
   int i, j;

   Q_data_TYP *data;



   data = (Q_data_TYP *)calloc(1, sizeof(Q_data_TYP));
   data->G = copy_bravais(G);
   data->pres = copy_mat(pres);
   data->l_option = l_option;

   /* calculate the Z-classes */
   data->INZ = (QtoZ_TYP *)calloc(1, sizeof(QtoZ_TYP));
   data->Z = q2z(G, &data->Z_no, 0, data->INZ, 1, TRUE);
   cleanup_prime();

   /* affine classes */
   data->all = 0;
   data->aff_no = (int *)calloc(data->Z_no, sizeof(int));
   data->WORDS = (int ****)calloc(data->Z_no, sizeof(int ***));
   data->NUMBER_OF_WORDS = (int **)calloc(data->Z_no, sizeof(int *));
   data->N = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   data->X = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   data->names = (MP_INT **)calloc(data->Z_no, sizeof(MP_INT *));
   if (!data->l_option)
      data->names_int = (int **)calloc(data->Z_no, sizeof(int *));
   data->gen_inv = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   data->coz = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   data->norm_inv = (matrix_TYP ***)calloc(data->Z_no, sizeof(matrix_TYP **));
   data->X_2_inv = (matrix_TYP **)calloc(data->Z_no, sizeof(matrix_TYP *));
   data->coho_size = (int *)calloc(data->Z_no, sizeof(int));
   data->stab_coz = (matrix_TYP ****)calloc(data->Z_no, sizeof(matrix_TYP ***));
   data->stab_gen_no = (int **)calloc(data->Z_no, sizeof(int *));
   data->first_aff = (int *)calloc(data->Z_no, sizeof(int));
   data->list_of_names = (int **)calloc(data->Z_no, sizeof(int *));
   for (i = 0; i < data->Z_no; i++){
      cen_to_norm(data->Z[i]);
      data->first_aff[i] = data->all;
      data->WORDS[i] = (int ***)calloc(MIN_SPEICHER, sizeof(int **));
      data->NUMBER_OF_WORDS[i] = (int *)calloc(MIN_SPEICHER, sizeof(int));
      data->gen_inv[i] = (matrix_TYP **)calloc(data->Z[i]->gen_no, sizeof(matrix_TYP *));
      data->coz[i] = all_cocycles(pres, data->Z[i], &data->aff_no[i], data->gen_inv[i],
                                  &data->X[i], &data->names[i], &data->WORDS[i],
                                  &data->NUMBER_OF_WORDS[i], &data->N[i], &data->coho_size[i],
                                  &data->list_of_names[i], data->l_option);
      data->all += data->aff_no[i];
      data->norm_inv[i] = (matrix_TYP **)calloc(data->Z[i]->normal_no, sizeof(matrix_TYP *));
      for (j = 0; j < data->Z[i]->normal_no; j++){
         data->norm_inv[i][j] = mat_inv(data->Z[i]->normal[j]);
      }
      data->X_2_inv[i] = mat_inv(data->X[i][2]);
      data->stab_coz[i] = (matrix_TYP ***)calloc(data->aff_no[i], sizeof(matrix_TYP **));
      data->stab_gen_no[i] = (int *)calloc(data->aff_no[i], sizeof(int));
      for (j = 0; j < data->aff_no[i]; j++){
         data->stab_coz[i][j] = stab_coz(data->WORDS[i][j], data->NUMBER_OF_WORDS[i][j],
                                         data->Z[i]->normal, data->norm_inv[i],
                                         data->Z[i]->normal_no, j, &data->stab_gen_no[i][j]);
      }
   }

   /* put cocycles and pointgroup together */
   data->aff = (bravais_TYP ***)calloc(data->Z_no, sizeof(bravais_TYP **));
   for (i = 0; i < data->Z_no; i++){
      data->aff[i] = (bravais_TYP **)calloc(data->aff_no[i], sizeof(bravais_TYP *));
      if (!data->l_option)
         data->names_int[i] = (int *)calloc(data->aff_no[i], sizeof(int));
      for (j = 0; j < data->aff_no[i]; j++){
         if (!data->l_option)
            data->names_int[i][j] = mpz_get_ui(&data->names[i][j]);
         data->aff[i][j] = extract_r(data->Z[i], data->coz[i][j]);
      }
   }
   return(data);
}




/* ---------------------------------------------------------------------------- */
/* Write the number of Z-classes and affine classes in each Z-class to stdout!  */
/* If printflag, we write the Z-classes and affine classes in files.            */
/* data->aff, data->aff_no, data->Z, data->Z_no have to be filled correctly     */
/* and completely                                                               */
/* ---------------------------------------------------------------------------- */
void put_Q_data(Q_data_TYP *data,
                char *groupname,
                int printflag)
{
   int i, j;

   char  comment[1000],
 	 file[1000];

   matrix_TYP *form, *id;

   id = init_mat(data->G->dim, data->G->dim, "1");

   if (printflag){
      for (i = 0; i < data->Z_no; i++){
         sprintf(comment, "%d-th Z-class of %s", i+1, groupname);
         sprintf(file, "%s.%d",groupname, i+1);
         put_bravais(data->Z[i], file, comment);

         form = rform(data->Z[i]->gen, data->Z[i]->gen_no, id, 100);
         sprintf(comment, "form of the %d-th Z-class of %s", i+1, groupname);
         sprintf(file, "form.%s.%d",groupname, i+1);
         put_mat(form, file, comment,0);
         free_mat(form);

         for (j = 0; j < data->aff_no[i]; j++){
            sprintf(comment, "%d-th affine class of the %d-th Z-class of %s",
		    j+1, i+1, groupname);
            sprintf(file, "%s.%d.%d", groupname, i+1, j+1);
            put_bravais(data->aff[i][j], file, comment);
         }
      }
   }

   printf("There are %i Z-Classes with ", data->Z_no);
   for (i = 0; i < data->Z_no; i++){
      printf("%i ",data->aff_no[i]);
   }
   printf("Space Groups!\n");

   free_mat(id);
}



/* ---------------------------------------------------------------------------- */
/* Free the structure Q_data_TYP                                                */
/* data has to be filled correctly and completely                               */
/* ---------------------------------------------------------------------------- */
void free_Q_data(Q_data_TYP *data)
{
   int i, j, k, Nanz;


   if (data->first_aff){
      free(data->first_aff);
   }
   if (data->names_int){
      for (i = 0; i < data->Z_no; i++){
         if (data->names_int[i] != NULL)
            free(data->names_int[i]);
      }
      free(data->names_int);
   }
   if (data->list_of_names){
      for (i = 0; i < data->Z_no; i++){
         if (data->list_of_names[i] != NULL)
            free(data->list_of_names[i]);
      }
      free(data->list_of_names);
   }
   if (data->stab_coz){
      for (i = 0; i < data->Z_no; i++){
         for (j = 0; j < data->aff_no[i]; j++){
            for (k = 0; k < data->stab_gen_no[i][j]; k++){
               free_mat(data->stab_coz[i][j][k]);
            }
            free(data->stab_coz[i][j]);
         }
         free(data->stab_coz[i]);
      }
      free(data->stab_coz);
   }
   if (data->stab_gen_no){
      for (i = 0; i < data->Z_no; i++){
         if (data->stab_gen_no[i] != NULL)
            free(data->stab_gen_no[i]);
      }
      free(data->stab_gen_no);
   }
   if (data->norm_inv){
      for (i = 0; i < data->Z_no; i++){
         for (j = 0; j < data->Z[i]->normal_no; j++)
            free_mat(data->norm_inv[i][j]);
         free(data->norm_inv[i]);
      }
      free(data->norm_inv);
   }
   if (data->gen_inv){
      for (i = 0; i < data->Z_no; i++){
         for (j = 0; j < data->Z[i]->gen_no; j++)
            if (data->gen_inv[i][j])
               free_mat(data->gen_inv[i][j]);
         free(data->gen_inv[i]);
      }
      free(data->gen_inv);
   }
   if (data->N){
      for (i = 0; i < data->Z_no; i++){
         if (data->N[i] != NULL){
            Nanz = data->Z[i]->cen_no + data->Z[i]->normal_no;
            for (j = 0; j < Nanz; j++)
               free_mat(data->N[i][j]);
            free(data->N[i]);
         }
      }
      free(data->N);
   }
   if (data->WORDS && data->NUMBER_OF_WORDS && data->X){
      for (i = 0; i < data->Z_no; i++){
         if (data->X[i][0]->cols >= 1){
            for (j = 0; j < data->aff_no[i]; j++){
               for (k = 0; k < data->NUMBER_OF_WORDS[i][j]; k++){
                  free(data->WORDS[i][j][k]);
               }
               free(data->WORDS[i][j]);
            }
         }
         free(data->WORDS[i]);
         free(data->NUMBER_OF_WORDS[i]);
      }
      free(data->WORDS);
      free(data->NUMBER_OF_WORDS);
   }
   if (data->X){
      for (i = 0; i < data->Z_no; i++){
         for (j = 0; j < 3; j++){
            if (data->X[i][j] != NULL) free_mat(data->X[i][j]);
         }
         free(data->X[i]);
      }
      free(data->X);
   }
   if (data->X_2_inv){
      for (i = 0; i < data->Z_no; i++){
         if (data->X_2_inv[i] != NULL) free_mat(data->X_2_inv[i]);
      }
      free(data->X_2_inv);
   }
   if (data->names)
     {
     for (i = 0; i < data->Z_no; i++)
       {
       if (data->names[i] != NULL)
          {
          for (j = 0; j < data->aff_no[i]; j++)
             mpz_clear(data->names[i] + j);
          free(data->names[i]);
          }
       }
     free(data->names);
     }
   if (data->coho_size){
      free(data->coho_size);
   }
   free_QtoZ(data->INZ, 0);
   for (i = 0; i < data->Z_no; i++){
      for (j = 0; j < data->aff_no[i]; j++){
         free_bravais(data->aff[i][j]);
         free_mat(data->coz[i][j]);
      }
      free_bravais(data->Z[i]);
      free(data->aff[i]);
      free(data->coz[i]);
   }
   free(data->aff_no);
   free(data->aff);
   free(data->coz);
   free(data->Z);
   free_bravais(data->G);
   free_mat(data->pres);
   free(data);
}










