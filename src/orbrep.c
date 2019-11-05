#include"typedef.h"


main (int argc, char *argv[])
{

	matrix_TYP **Mat, *vecs, *represent;
        bravais_TYP *G;
        int orb_anz;
        int *subdiv;
        int *c, *f;
        int i,j, anz;

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
           printf("usage:   orbrep 'file1' 'file2',\n");
           printf("where  'file1' contains a matrix_TYP and\n");
           printf("       'file2' contains a bravais_TYP\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = mget_mat (FILENAMES[0], &anz);
        vecs = Mat[0];
        G = get_bravais(FILENAMES[1]);
        subdiv = orbit_subdivision(vecs, G, &orb_anz);
        represent = init_mat(orb_anz, G->gen[0]->cols, "");
        c = (int *)malloc(orb_anz *sizeof(int));
        f = (int *)malloc(orb_anz *sizeof(int));
        for(i=0;i<orb_anz;i++)
           c[i] = 0;
        for(i=0;i<vecs->rows;i++)
        {
          j = subdiv[i]-1;
          if(c[j] == 0)
            f[j] = i;
          c[j]++;
        }
        for(i=0;i<orb_anz;i++)
        {
           for(j=0;j<G->gen[0]->cols;j++)
             represent->array.SZ[i][j] = vecs->array.SZ[f[i]][j];
        }
        put_mat(represent, NULL, "representatives of the orbits", 0);
        printf(" The length of the i-th orbit is:\n");
        for(i=0;i<orb_anz;i++)
          printf("%d  ", 2  *c[i]);
        printf("\n");


   exit(0);
}
