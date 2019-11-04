#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <sort.h>
#include <bravais.h>
#include <base.h>
#include <longtools.h>
#include <datei.h>
#include <presentation.h>


#define DEBUG FALSE


static matrix_TYP *search_inverse(bahn *s,
                                  matrix_TYP *x)
{

  int i;

  for (i=0;i<s->length;i++)
     if (s->representatives[i] == x)
        return s->rep_invs[i];

   return (matrix_TYP *) -1;

} /* search_inverse(...) */



static void con_reduce(int *w){

  int i = 1;

  while (i){
    if (w[1] == -w[w[0]]){
      for (i=1;i<w[0];i++) w[i] = w[i+1];
      w[0] -= 2;
    }
    else{
       i = 0;
    }
  }

  return;

} /* con_reduce */


void normalize_word(int *w){

  int i,
      j;

  for (i=w[0];i>1;i--){
    if (w[i] == -w[i-1]){
      for (j=i+1;j<=w[0];j++) w[j-2] = w[j];
      w[0] -= 2;
      i = w[0] + 1;
    }
  }

  return;

} /* normalize_word(int *w) */



static void order(int **w,
                  int *l)
{

  int i,
     *P,
      good = 1;

  while (good){
    good = 0;

    for (i=0;i<l[0]-1;i++){
      if (w[i][0] > w[i+1][0]){
	P = w[i];
        w[i] = w[i+1];
        w[i+1] = P;
        good = 1;
      }
    }
  }


  return;

} /* order(...) */

static int skip_sub_relator(int *w1,int *w2){

  int i,j,n,*wn;

  n = w1[0]/2 + 1;

  for (i=0;i<(w2[0] - n +1);i++){
    if (!memcmp(w1+1,w2+1+i,n*sizeof(int))){
      wn = (int *) malloc((1+w1[0]+w2[0]) * sizeof(int));
      wn[0] = w1[0]+w2[0];
      memcpy(wn+1,w2+1,i * sizeof(int));
      for (j=1;j<=w1[0];j++)
	wn[j+i] = -w1[w1[0] - j + 1];
      memcpy(wn+i+1+w1[0],w2+i+1,(w2[0] - i) * sizeof(int));

      normalize_word(wn);
      if (wn[0] >= w2[0]){
	fprintf(stderr,"error\n");
        exit(3);
      }

      memcpy(w2,wn,(wn[0] + 1) * sizeof(int));
      free(wn);
      return TRUE;
    }

  }

  return FALSE;
} /* skip_sub_relator(...) */


static int skip_sub_relator_inv(int *w1,int *w2){

  int i,
      j;

  static int wi[11];

  if (w1[0] > 10) return FALSE; 

  for (i=1,j=w1[0] ; i<=w1[0] ; i++,j--)
    wi[i] = -w1[j];
  wi[0] = w1[0];

  return skip_sub_relator(wi,w2);

} /* skip_sub_relator_inv(...) */



void simplify_presentation(int **w,
                           int *l){
  int i,
      j,
      good = 1;

  while (good){
    good = 0;

    /* search for trivial relators */
    for (i=0;i<l[0];i++){
      if (w[i][0] == 0){
	free(w[i]);
        l[0]--;
        w[i] = w[l[0]];
        good = 1;
        i = l[0] + 1;
      }
    }

    /* oder the relators according to their length */
    order(w,l);


    /* see if one of the short relators is a subrelator of a longer one */
    for (i=0;i<l[0] - 1; i++){
      for (j=l[0]-1;j>i ; j--){
	if (w[i][0] > 0 &&
            w[j][0] > 0){
          if(skip_sub_relator(w[i],w[j]) || skip_sub_relator_inv(w[i],w[j]))
	    good = TRUE;
        }
      }
    }

    for (i=0;i<l[0];i++){
      normalize_word(w[i]);
      con_reduce(w[i]);
    }

    j = 0 ;
    for (i=0;i<l[0];j+=w[i][0] , i++);

    if (DEBUG)
      fprintf(stderr,"%d relators of total length %d\n",i,j);

  }


  return;

} /* simplify_presentation(...) */


static int *relator3(int k,
                     int *w1,
                     int *w2){

  int i,
      j,
     *res;


  res = (int *) calloc( w1[0] + w2[0] + 2, sizeof(int));
  res[0] = w1[0] + w2[0] + 1;

  memcpy(res+2,w1+1,w1[0] * sizeof(int));
  res[1] = k+1;

  for (i=1,j=res[0];i<=w2[0];j--,i++)
    res[j] = -w2[i];

  return res;

} /* relator3 (....) */

static int *relator2(int *w1,
                     int *wz,
                     int *w2){

  int i,
      j,
     *res;


  res = (int *) calloc( w1[0] + w2[0] + wz[0] + 1 , sizeof(int));
  res[0] = w1[0] + w2[0] + wz[0];

  memcpy(res+1,w1+1,w1[0] * sizeof(int));
  memcpy(res+1+w1[0],wz+1,wz[0] * sizeof(int));

  for (i=1,j=res[0];i<=w2[0];j--,i++)
    res[j] = -w2[i];

  return res;

} /* relator2 (....) */


/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ void put_presentation(int **w,
@                       int l,
@                       bravais_TYP *G,
@                       char *O)
@
@
@
@
@--------------------------------------------------------------------------
@
***************************************************************************/
void put_presentation(int **w,
                      int l,
                      bravais_TYP *G,
                      const char *O){

  int i;


  if (strchr(O,'G')){
    if (strchr(O,'4')){
      printf("F := FreeGroup(%d);\n",G->gen_no);
      printf("g := GeneratorsOfGroup(F);");
    }

    printf("rel := [\n");
  }

  for (i=0;i<l-1;i++){
    put_word(w[i],O);

    if (strchr(O,'G')){
      printf(",\n");
    }
  }

  put_word(w[i],O);
  if (strchr(O,'G')){
    printf("];\n");
  }


  if (strchr(O,'G')){

    printf("G := F/ rel;\n");

    if (strchr(O,'V')){
       printf("p := PresentationFpGroup(G);\n");
       printf("TzGoGo(p);\n");
       printf("G := FpGroupPresentation(p);\n");
    }

    if (strchr(O,'S')){
      printf("Size(G);\n");
    }

    if (strchr(O,'Q')){
      printf("quit;\n");
    }
  }

  return;

}  /* put_presentation(...) */




static int *get_word_to_generator(bahn *A,
                                  matrix_TYP *h)
{

  int i;

  for (i=0;i<A->length && mat_comp(A->representatives[i],h) ; i++);

  if (i == A->length){
    fprintf(stderr,"error\n");
    exit(3);
  }

  return A->words[i];

} /* get_word_to_generator */



int *mul_word(int *w1,
              int *w2){

   int *w;


   w = (int *) malloc((w1[0]+w2[0]+1) * sizeof(int));

   w[0] = w1[0]+w2[0];

   memcpy(w+1,w1+1,w1[0]*sizeof(int));
   memcpy(w+1+w1[0],w2+1,w2[0]*sizeof(int));

   return w;

}


static int needed(matrix_TYP **VL,
                  matrix_TYP *x,
                  matrix_TYP **GEN){

  int i = 0,
      j,
      k,
      RES = FALSE;

  matrix_TYP *y;


  i=0;
  while(VL[i]){
     y = mat_mul(x,VL[i]);

     j=0;
     while(VL[j] && y){
       if (!mat_comp(VL[j],y)){
          free_mat(y);
          y = NULL;
       }
       j++;
     }
     if (y){
        VL[j] = y;
        RES = TRUE;
     }

     i++;
  }

  if (RES){

     i=0;
     while(GEN[i]) i++;
     GEN[i] = x;

     i = 0;
     while (VL[i]){

       k = 0;
       while(GEN[k]){
         y = mat_mul(GEN[k],VL[i]);

         j=0;
         while(VL[j] && y){
           if (!mat_comp(VL[j],y)){
             free_mat(y);
             y = NULL;
	   }
           j++;
         }
         if (y){
           VL[j] = y;
         }
         k++;
       }
       i++;
     }
  }

  return RES;

} /* needed (.....) */


/*************************************************************************
@
@-------------------------------------------------------------------------
@
@
@-------------------------------------------------------------------------
@
**************************************************************************/
static void con_relations(bravais_TYP *G,
                          matrix_TYP **GGENINV,
                          matrix_TYP **GEN,
                          matrix_TYP **GENINV,
                          int **GW,
                          bahn **s,
                          int level,
                          int ***w,
                          int *l,
                          int *speicher,
                          matrix_TYP **GENU,
                          int **GENUW,
                          int GENU_NO)
{


   matrix_TYP **VS,
              **GENUH,
               *y,
               *x;


   int i,
       j,
       k,
      *w2 = NULL,
      *w3,
     **GENUWH,
       found = GENU_NO;



   VS = (matrix_TYP **) calloc(s[level]->length+1,sizeof(matrix_TYP *));
   VS[0] = s[level]->orbit[0];

   GENUWH = (int **) malloc((GENU_NO+s[level]->length) * sizeof(int *));
   GENUH = (matrix_TYP **) calloc((GENU_NO+s[level]->length),
                                   sizeof(matrix_TYP*));
   memcpy(GENUH,GENU,GENU_NO * sizeof(matrix_TYP *));
   memcpy(GENUWH,GENUW,GENU_NO * sizeof(int *));

   /* see which conjugate we will need */
   for (i=0;i<found;i++){
     k = TRUE;
     for (j=0;j<s[level]->gen_no && k;j++){

        x = mat_kon(GEN[j],GENUH[i],GENINV[j]);

        if (needed(VS,x,GENUH+GENU_NO)){
           GENUWH[found] = relator2(GW[j],GENUWH[i],GW[j]);

           if (DEBUG){
             y = mapped_word(GENUWH[found],G->gen,GGENINV);

             if (mat_comp(x,y)){
	       fprintf(stderr,"error in pres\n");
               exit(3);
	     }
             free_mat(y);
           }

           found++;
        }
        else{
           free_mat(x);
        }

        for (k=0;VS[k];k++);
        k = (k==s[level]->length);

     }
   }


   if (DEBUG)
     fprintf(stderr,"GENU_NO %d found %d level %d\n",GENU_NO,found,level);

   /* conjugate relations */
   for (i=0;i<found;i++){
     for (j=0;j<s[level]->gen_no;j++){

       x = mat_kon(GEN[j],GENUH[i],GENINV[j]);

       /* find a word in the generators, and check a bit */
       if (!is_element(x,G,s,&w2)){
          fprintf(stderr,"error in pres\n");
	  exit(3);
       }

       if (DEBUG){
          y = mapped_word(w2,G->gen,GGENINV);

          if (mat_comp(x,y)){
	     fprintf(stderr,"error in pres\n");
             exit(3);
	  }
          free_mat(y);
       }

       w3 = relator2(GW[j],GENUWH[i],GW[j]);

       if (DEBUG){
          y = mapped_word(w3,G->gen,GGENINV);

          if (mat_comp(x,y)){
	     fprintf(stderr,"error in pres\n");
             exit(3);
	  }
          free_mat(y);
       }

       if ( *l == *speicher){
         *speicher += 256;
         *w = (int **) realloc(*w,*speicher * sizeof(int *));
       }


       k = 0;
       w[0][*l] = relator2(&k,w3,w2);

       normalize_word(w[0][*l]);
       con_reduce(w[0][*l]);

       if (DEBUG){
         y = mapped_word(w[0][*l],G->gen,GGENINV);
         Check_mat(y);
	 if (!y->flags.Diagonal){
	   fprintf(stderr,"error in pres\n");
           exit(3);
	 }
         free_mat(y);
       }

       if (w[0][*l][0] == 0){
	 free(w[0][*l]);
       }
       else{
         /* put_word(w[0][*l],"M"); */
         (*l)++;
       }

       free(w3);
       free(w2); w2 = NULL;
       free_mat(x);
     }
   }


   for (i=GENU_NO;i<found;i++){
      free_mat(GENUH[i]);
      free(GENUWH[i]);
   }
   free(GENUH);
   free(GENUWH);

   for (i=1;i<s[level]->length;i++)
     if (VS[i]) free_mat(VS[i]);
   free(VS);

   return;

}



matrix_TYP *pres_on_many_generators(bravais_TYP *G,
				    int *OPT){


    bravais_TYP *H;
    bahn **s;
    matrix_TYP **basis;
    int i,j;
    int k,l;
    int zeile;
    matrix_TYP *RES;
    int MAX = G->gen_no + 2;
    int *w;

    s = NULL;
    H = copy_bravais(G);

    basis = get_base(H);


    if (H->gen_no > H->dim){
	H->order = red_gen(H,basis,&s,0);
	for (i=0;i<H->dim;i++){
	    free_bahn(s[i]);
	    free(s[i]);
	}
	free(s);
    }

    s = strong_generators(basis,H,TRUE);

    /* the situation is: H = G as a group, and the generators of H form a subset of the generators of G */

    /* calculate a presentation for H on the generators of H */
    RES = pres(s,H,OPT);

    zeile = RES->rows;
    real_mat(RES,RES->rows + G->gen_no  , RES->cols);

    /* add relations of the form g_i = word(h_j) for those g_i which are not in the generators of H, */
    /* but represent the g_i by i+MAX */
    /* mathemetically this is an tietze tranformation of the third (?) kind, add a generator */
    /* <h_i|rs> \cong <h_i,g| rs,g=w(h_i) > */
    for (i=0;i<G->gen_no;i++){

	/* is_element calculates a word such that g_i = w(h_j) */
	w = NULL;
	if(!is_element(G->gen[i],H,s,&w)){
	    fprintf(stderr,"error in pres_on_many_generators(...)\n");
	    exit(3);
	}

	if (RES->cols <= w[0]){
	    real_mat(RES,RES->rows,w[0]+1);
	}

	/* the relation is now g_i^(-1) * w(h_J) = 1 */
	memcpy(RES->array.SZ[zeile],w,(w[0]+1) * sizeof(int));
	RES->array.SZ[zeile][0] = -i - 1 - MAX;

	/* for sanitary reasons remove those relations which look like g_i * g_i^(-1) */
	if (w[0] == 1 && w[1] == i+1){
	    RES->array.SZ[zeile][0] = 0;
	    RES->array.SZ[zeile][1] = 0;
	}
	else{
	    zeile++;
	}

	free(w);
    }

    if (zeile < RES->rows)
	real_mat(RES,zeile,RES->cols);


    /* translate this presentation in the generators of G */
    /* on strange thing : to be able to do this in one go, we represent the j-th generator of G */
    /* by the number j+MAX temporaly */
    for (i=0; i< H->gen_no ; i++){

	/* search for a j with H->gen[i] == G->gen[j] */
	for (j=0; j< G->gen_no && mat_comp(H->gen[i],G->gen[j]) ; j++);

	/* replace every occurence of i+1 or -(i+1) in RES with (j+1)+MAX or -(j+1)-MAX respectively */
	for (k=0;k<RES->rows;k++){
	    for (l=0;l<RES->cols;l++){
		if (RES->array.SZ[k][l] == i+1)
		    RES->array.SZ[k][l] = j + 1 + MAX;
		else if (RES->array.SZ[k][l] == -i-1)
		    RES->array.SZ[k][l] = -j - 1 - MAX;
	    }
	}
    }

    /* now RES contains relations for G with the strange thing that the number j generator has been
       replaced with the generator j+MAX temporaly if it is also a generator of H */
    for (k=0;k<RES->rows;k++){
	for (l=0;l<RES->cols;l++){
	    if (RES->array.SZ[k][l] > MAX)
		RES->array.SZ[k][l] -= MAX;
	    else if (RES->array.SZ[k][l] < -MAX)
		RES->array.SZ[k][l] += MAX;
	}
    }



    /* clean up */
    for (i=0;i<H->dim;i++){
	free_mat(basis[i]);
	free_bahn(s[i]);
	free(s[i]);
    }
    free(basis);
    free(s);
    free_bravais(H);

    return RES;


} /* pres_on_many_generators(...) */


/************************************************************************
@
@------------------------------------------------------------------------
@
@ matrix_TYP *pres(bahn **s,
@                  bravais_TYP *G,
@                  int *OPT)
@
@
@
@------------------------------------------------------------------------
@
*************************************************************************/
matrix_TYP *pres(bahn **s,
                 bravais_TYP *G,
                 int *OPT){

  matrix_TYP *M,
             *x,
             *y,
             *g,
            **GENU,
            **GENUINV,
            **GENINV = NULL;

  static int has_been_called_recursively;

  int i,
      ii,
      j,
      k,
      l = 0,
      GENU_NO = 0,
    **GENUW,
    **GW,
      speicher = 256 - 8,
     *w2 = NULL,
    **w;

  if ((!has_been_called_recursively && G->gen_no > G->dim) || s == NULL){
      /* many gerators, treat it differently */
      /* and prevent using users the form pres_on_many_generators just because they are too lazy to calculate */
      /* s, and calculate it for them */
      has_been_called_recursively = TRUE;
      return pres_on_many_generators(G,OPT);
  }

  has_been_called_recursively = FALSE;

  GENU = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
  GENUINV = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
  GENUW = (int **) malloc(sizeof(int *));

  GENINV = (matrix_TYP **) calloc(G->gen_no,sizeof(matrix_TYP*));

  for (i=0;i<G->gen_no;i++)
    GENINV[i] = long_mat_inv(G->gen[i]);

  /* alloc space for the result */
  w = (int **) calloc(speicher , sizeof(int*));

  s[0]->gen_no = G->gen_no;

  for (i=G->dim-1;i>=0;i--){
    for (j=0;j<s[i]->length;j++){

      /* relations of the first kind */
      if (i==0){
          ii=G->gen_no;
      }
      else{
          ii = G->gen_no + s[i]->gen_no;
      }
      for (k=0;k<ii;k++){

        if ( l == speicher){
           speicher += 256;
           w = (int **) realloc(w,speicher * sizeof(int*));
        }

        if (k < G->gen_no){
           g = G->gen[k];
        }
        else{
           g = s[i]->generators[k - G->gen_no];
        }

        x = mat_mul(g,s[i]->representatives[j]);

        if (!is_element(x,G,s,&w2)){
	  fprintf(stderr,"error in pres\n");
	  exit(3);
	}

        if (DEBUG){
          y = mapped_word(w2,G->gen,GENINV);

          if (mat_comp(x,y)){
	    fprintf(stderr,"error in pres\n");
            exit(3);
	  }
          free_mat(y);

        }


        if (k<G->gen_no){
          w[l] = relator3(k,s[i]->words[j],w2);
        }
        else{
          w[l] = relator2(get_word_to_generator(s[i],g),s[i]->words[j],w2);
	}

        normalize_word(w[l]);
        con_reduce(w[l]);

        if (DEBUG){
          y = mapped_word(w[l],G->gen,GENINV);
          Check_mat(y);
          if (!y->flags.Diagonal){
	    fprintf(stderr,"error in pres\n");
            exit(3);
	  }
          free_mat(y);
	}

        if (w[l][0] == 0){
	  free(w[l]);
	}
	else{
          /* put_word(w[l],"M"); */
          l++;
        }

        free_mat(x);
        free(w2); w2 = NULL;

      }


    }

      if (i<G->dim -1 ){
        /* relations of the second kind, conjugation */
        GENU = (matrix_TYP **) realloc(GENU,(GENU_NO+1+s[i+1]->gen_no)
                                            *sizeof(matrix_TYP *));
        GENUW = (int **) realloc(GENUW,(GENU_NO+1+s[i+1]->gen_no)
                                            *sizeof(int*));

        for (k=0;k<s[i+1]->gen_no;k++){
           GENU[GENU_NO + k] = s[i+1]->generators[k];
           GENUW[GENU_NO + k] = get_word_to_generator(s[i+1],
                                s[i+1]->generators[k]);
        }
        GENU_NO += s[i+1]->gen_no;

        GW = (int **) malloc(s[i]->gen_no * sizeof(int*));
        if (i==0){

           for (k=0;k<G->gen_no;k++){
             GW[k] = (int *) malloc(2*sizeof(int));
             GW[k][0] = 1;
             GW[k][1] = k+1;
           }

           con_relations(G,GENINV,G->gen,GENINV,
                         GW,s,i,
                         &w,&l,&speicher,
                         GENU,GENUW,GENU_NO);

           for (k=0;k<G->gen_no;k++)
             free(GW[k]);
        }
        else{

           GENUINV = (matrix_TYP **) realloc(GENUINV,
                                    s[i]->gen_no*sizeof(matrix_TYP*));
           for (k=0;k<s[i]->gen_no;k++){
             GW[k] = (int *) get_word_to_generator(s[i],s[i]->generators[k]);
             GENUINV[k] = search_inverse(s[i],s[i]->generators[k]);
           }

           con_relations(G,GENINV,s[i]->generators,GENUINV,
                         GW,s,i,
                         &w,&l,&speicher,
                         GENU,GENUW,GENU_NO);
        }

        free(GW);

      }


  }

  /* put_presentation(w,l,G,"G4"); */

  simplify_presentation(w,&l);

  if (DEBUG){
    for (i=0;i<l;i++){
      y = mapped_word(w[0],G->gen,GENINV);
      Check_mat(y);
      if (!y->flags.Diagonal){
	fprintf(stderr,"error in pres\n");
	exit(3);
      }
      free_mat(y);
    }
  }

  if (GENINV){
    for (i=0;i<G->gen_no;i++)
      if (GENINV[i]) free_mat(GENINV[i]);
    free(GENINV);
  }


  if (OPT && OPT[0])
     put_presentation(w,l,G,"GSV4");

  j=0;
  for (i=0;i<l;i++) if (j<w[i][0]) j = w[i][0];
  M = init_mat(l,j,"i");
  for (i=0;i<l;i++){
    memcpy(M->array.SZ[i],w[i]+1,w[i][0]*sizeof(int));
    free(w[i]);
  }
  free(w);
  free(GENU);
  free(GENUINV);
  free(GENUW);

  /* put_mat(M,0,0,0); */

  return M;
} /* pres(....) */


