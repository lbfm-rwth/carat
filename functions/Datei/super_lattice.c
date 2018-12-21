#include "typedef.h"
#include "datei.h"
#include "matrix.h"
#include "getput.h"

lattice_element **super_lattice(char *symb,int dim,int almost,int zclass,
		int *no, int OPTION)
{

  int i, j, pos, found, c;

   FILE *infile;

   lattice_element **RES,
                    *TMP;

   /* saves space on the stack */
   static char filename[1024], format[1024];

   /* get the appropiate filename */
   get_data_dir(format, "tables/lattices/dim%d/%s%s_%d_%d");
   sprintf(filename, format, dim, "reverse_", symb, almost, zclass);

   infile = fopen(filename,"r");

   /* save cores */
   if (infile == NULL){
     fprintf(stderr,"Couldn't find my input file %s\n",filename);
     exit(4);
   }

   c=fscanf(infile,"%d\n",no);

   RES = (lattice_element **) malloc(no[0] * sizeof(lattice_element *));

   for (i=0;i<no[0];i++){
      RES[i] = init_lattice_element();

      c=fscanf(infile,"%s %d %d\n",
      RES[i]->symbol,&RES[i]->almost,&RES[i]->zclass);
   }

   fclose(infile);

   if (OPTION){
      for (i=0;i<no[0];i++){

         found = FALSE;
         /* get the appropiate filename */
         sprintf(filename, format, dim, "lattice_",
		 RES[i]->symbol, RES[i]->almost, RES[i]->zclass);

         infile = fopen(filename,"r");
         c=fscanf(infile,"#%d\n",&pos);
         j = 0;
         while (!found && j<pos){
            TMP = fget_lattice_element(infile,FALSE);

            found = (strcmp(TMP->symbol,symb) == 0);
            found = found && (TMP->almost == almost) && (TMP->zclass == zclass);

            if (found){
               RES[i]->TR = TMP->TR;
               TMP->TR = NULL;
               RES[i]->N_orbits = TMP->N_orbits;
               TMP->N_orbits = 0;
            }

            j++;
            free_lattice_element(TMP);
         }
         fclose(infile);

         if (!found){
            fprintf(stderr,"ERROR: super_lattice\n");
            exit(3);
         }

         /* get the apropriate group */
         RES[i]->grp = brav_from_datei(RES[i]->symbol,
                                       RES[i]->almost,RES[i]->zclass);
      }
   }

   return RES;
}

