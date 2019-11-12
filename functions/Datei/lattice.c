#include "typedef.h"
#include "datei.h"


lattice_element **lattice(const char *symb,int dim,int almost,int zclass,int *no,
                          int OPTION)
{
   lattice_element **res;

   FILE *infile;

   int i;

   /* saves space on the stack */
   static char filename[1024];

   /* get the appropiate filename */
   sprintf(filename, "%s/lattices/dim%d/lattice_%s_%d_%d", get_data_dir(), dim, symb, almost, zclass);

   infile = fopen(filename,"r");

   /* save cores */
   if (infile == NULL){
     fprintf(stderr,"Couldn't find my input file %s\n",filename);
     exit(4);
   }

   fscanf(infile,"#%d\n",no);

   if (INFO_LEVEL & 1){
      fprintf(stderr,"%d\n",no[0]);
   }

   /* now read the data for lattice */
   res = (lattice_element **) malloc(no[0] * sizeof(lattice_element *));

   /* read the subgroups */
   for (i=0;i<no[0];i++){

      /* get the element of lattice */
      res[i] = fget_lattice_element(infile,OPTION);

      /* output for debugging */
      if (INFO_LEVEL & 1){
         fput_lattice_element(res[i],NULL);
      }
   }

   fclose(infile);

   return res;
}


