#include"typedef.h"
#include"bravais.h"
#include"base.h"

#include"voronoi.h"
#include"matrix.h"
#include"getput.h"
#include"symm.h"
#include"autgrp.h"
#include"reduction.h"

int INFO_LEVEL;
extern int SFLAG;

main (int argc, char *argv[])
{
   bravais_TYP *G,      /* the bravais group of H */
               *Gtr,
               *H;

   matrix_TYP *A,
              *tmp,
              *tmp2;

   voronoi_TYP **V;

   int Vanz,
       prime,
       i;


   read_header(argc, argv);

   if((FILEANZ < 1) || (is_option('h') && optionnumber('h')==0))
   {
     printf("Usage: Normalizer 'file1' ['file2' ['file3']] [-p=prime] [-o]\n");
     printf("\n");
     printf("file1: bravais_TYP containing the group G.\n");
     printf("file2: (OPTIONAL) bravais_TYP containing the transposed of G. (cf. Tr_bravais)\n");
     printf("file3: (OPTIONAL) matrix_TYP of a G-perfect form. (cf. First_perfect)\n");
     printf("\n");
     printf("Calculates a set of matrices which together with G generate the\n");
     printf("normalizer N_GL_n(Z) (G) of G in GL_n(Z).\n");
     printf("NOTE: the output echoes any input information about G, except input about\n");
     printf("generators of the normalizer.\n");
     printf("NOTE: The dimension of the space of invariant forms is a measure for the\n");
     printf("complexity of the algorithm. Up to degree 6 the only infeasible case are\n");
     printf("<I_6> and <-I_6>. Here the generators of the normalizer can be taken\n");
     printf("from `Bravais_cat' with family 1,1,1,1,1,1.\n");
     printf("\n");
     printf("Options:\n");
     printf("-b      :  The normalizer of the bravais group B(G) is calculated. With this\n");
     printf("           option the program is much faster. (The normalizer of G is a\n");
     printf("           subgroup of N_GL_n(Z) (B(G)). )\n");
     printf("-p=prime:  The determinants of the perfect forms are\n");
     printf("           calculated module prime. The default is 1949.\n");
     printf("-o      :  The G-perfect forms are given as additional output.\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
   }

   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   H = get_bravais(FILENAMES[0]);

   /* throw away the normalizer of H */
   for (i=0;i<H->normal_no;i++){
      free_mat(H->normal[i]);
   }
   if (H->normal != NULL && H->normal_no > 0) free(H->normal);
   H->normal = NULL;
   H->normal_no = 0;

   if (is_option('b')){
      G = copy_bravais(H);
   }
   else{
      G = bravais_group(H);
   }

   /* let's see whether we already got the formspace */
   if (G->form == NULL){
      G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
   }

   /* read the transposed group if it is given */
   if (FILEANZ > 1){
      Gtr = get_bravais(FILENAMES[1]);
   }
   else{
      Gtr = tr_bravais(G,1,FALSE);
   }

   /* read an G-perfect form if it is given */
   if (FILEANZ > 2){
      A = get_mat(FILENAMES[2]);
   }
   else{
      /* firstly calculate an positive definite G-invariant form */
      tmp2 = init_mat(G->dim,G->dim,"1");
      tmp = rform(G->gen,G->gen_no,tmp2,101);
      free_mat(tmp2);

      /* now calculate the trace bifo */
      tmp2 = trace_bifo(G->form,Gtr->form,G->form_no);
      A = first_perfect(tmp,G,Gtr->form,tmp2,&Vanz);
      free_mat(tmp2);
      free_mat(tmp);
   }
   prime = optionnumber('p');
   if(prime == 0)
    prime = 1949;
   V = normalizer(A, G, Gtr, prime, &Vanz);

   /* now we got G and it's normalizer, so see what we can do with H */
   if (is_option('b')){
      /* very easy in this case */
      H->normal = G->normal;
      H->normal_no = G->normal_no;
      G->normal = NULL;
      G->normal_no = 0;

   }
   else{
      H->normal = normalizer_in_N(H,G,&H->normal_no,FALSE);
   }

   put_bravais(H, NULL, "group with complete normalizer");

   if(is_option('o')){
      for(i=0;i<Vanz;i++)
         put_voronoi(V[i]);
   }

   /* cleaning up the memory */
   for(i=0;i<Vanz;i++){
      clear_voronoi(V[i]);
      free(V[i]);
   }
   free(V);
   free_mat(A);
   free_bravais(H);
   free_bravais(G);
   free_bravais(Gtr);

   /* some diagnostic for memory leakage */
   if (INFO_LEVEL & 12){
      pointer_statistics(0,0);
   }
   exit(0);

}

