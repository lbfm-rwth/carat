
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: Pres.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include"typedef.h"
#include"matrix.h"
#include"polyeder.h"
#include"sort.h"
#include"presentation.h"
#include"getput.h"

static corner_TYP *generate_corners(polyeder_TYP *Pol, int *anz);
static int **bary_list(polyeder_TYP *Pol);
static int* generate_Liste_with_check(polyeder_TYP *Pol,int **Liste_inv);
static int search_mat(int test, int *Liste, int anz);
static int search_corner(int s1, int s2, corner_TYP *corners);
static int satisfies_map_extra( polyeder_TYP *Pol, int w, matrix_TYP *vec, matrix_TYP *M);
static int is_on_wall(corner_TYP corn,int *wallungl,polyeder_TYP *Pol,int d,int w,matrix_TYP *M);
static int *map(polyeder_TYP *Pol, int i, matrix_TYP *M);
static int g_wall(int s2, corner_TYP cork,int * Liste_inv, polyeder_TYP *Pol);
static int create_new_wall(int s2, int s1, int *Liste_inv, int d, polyeder_TYP *Pol, corner_TYP *corners, int anz);
static void append_relation(int S,int *lhsproduct,int lhs_no,int *Liste);
static int recall(int rel, int* erzlist, int gen_no);
static matrix_TYP *pres_word_to_mat(int *word, int length, polyeder_TYP *Pol, int *Liste);
static int calc_order(matrix_TYP *A, matrix_TYP *B);


/*  Fkt, die eine Liste der Kanten (Schnitte zweier Waende) erstellt.
    merke: No. der Waende die zu dieser Kante gehoeren
           Liste der Vertices auf dieser Kante
	   das Baryzentrum dieser Kante
*/
static corner_TYP *generate_corners(polyeder_TYP *Pol, int *anz)
{
  corner_TYP	*erg;
  int 		i,j,k,l,m;
  int		count;
  int 		iok,jok;
  int 		d;
  
  d = Pol->vert[0]->dim;
  
  m = (Pol->wall_no)*(Pol->wall_no-1)/2;
  erg = (corner_TYP *)malloc(m * sizeof(corner_TYP));
  count = 0;
  m=0;
  for(i=0; i<Pol->wall_no; i++)
    {
      for(j=0; j<i; j++)
	{
	  erg[count].w[0] = i;
	  erg[count].w[1] = j;
	  erg[count].v = (int*)calloc((Pol->vert_no),sizeof(int));
	  for(k=0; k<Pol->vert_no; k++)
	    {
	      jok = 0;
	      iok = 0;
	      for(l=0; l<Pol->vert[k]->wall_no; l++)
		{
		  jok = jok || (Pol->vert[k]->wall[l] == j);
		  iok = iok || (Pol->vert[k]->wall[l] == i);
		}
	      if(iok && jok)
		{
		  erg[count].v[m] = k;
		  m++;
		}
	    }/**for(k<Pol->vert_no)**/
	  if(m >= d-2)
	    {
	      erg[count].v = (int*)realloc(erg[count].v,m*sizeof(int));
	      erg[count].v_anz = m;
	      count ++;
	      m=0;
	    }
	  else{
	    m = 0;
	  }
	}/**for(j<i)**/
    }/**for(i<Pol->wall_no)**/
  *anz = count;
  return(erg);
}/**generate_corners()**/

/*berechne baryzentrum der Waende von Pol */
static int **bary_list(polyeder_TYP *Pol)
{
  int **erg;
  int *v_wand;	 	/*Nummern der Ecken der Wand*/
  int v_anz; 		/*Anzahl der Ecken der Wand*/
  int i, j, k, a, b;
  int *list;
  int kgv;
  int d;
  
  d = Pol->vert[0]->dim;
  erg = (int **)malloc(Pol->wall_no * sizeof(int*));
  for(i=0; i<Pol->wall_no; i++){
    erg[i] = (int *)malloc(d * sizeof(int));
    for(k=0; k<d; k++)
      erg[i][k] = 0;
  }
  
  for(k=0; k<Pol->wall_no; k++){
    
    /* bestimme v_wand und v_anz */
    v_anz = 0;
    v_wand = (int *) calloc(Pol->vert_no , sizeof(int));  
    
    /* bestimme Eckenmenge der Wand*/
    for(a=0; a<Pol->vert_no; a++)     
      { 
	for(b=0; b<Pol->vert[a]->wall_no; b++)
	  {
	    if(Pol->vert[a]->wall[b] == k)
	      {	
	      	v_wand[v_anz] = a;
              	v_anz++;
	      }
	  }
      }
    v_wand = (int*)realloc(v_wand,v_anz * sizeof(int));
	
    /*bestimme Baryzentrum der Wand*/
    list = (int*) calloc(v_anz ,sizeof(int));
    for(i=0; i<v_anz; i++) 		
      {
	/*v_wand[i]= Originalnummer der i-tn Ecke im Ausgangspolyeder*/
	list[i] = Pol->vert[v_wand[i]]->v[d-1];
      }
    kgv = KKGV(list, v_anz);	/*funktion KKGV extern */ 
    
      for(j=0; j<d; j++){
	for(i=0; i<v_anz; i++)
          {
	    erg[k][j] +=  (kgv/list[i])*(Pol->vert[v_wand[i]]->v[j]); 
    	    /* teile noch durch letzte Komponente, diese ist Anzahl*kgv */
          }
      }
      free(list);
  }
  free(v_wand);

  return(erg);
}/**bary_list()**/

/*
  Liste[i] = j -> Pol->wall[i] hat Seitentrafo = Pol->wall[j-1]->mat falls j>0 
  und Seitentrafo = (Pol->wall[j-1]->mat)^-1 ; j<0.
  */
static int* generate_Liste_with_check(polyeder_TYP *P,int **Liste_inv)
{
  int 			*Liste;
  int 			i,j,s,a;
  matrix_TYP 		*M_inv;
  int 			*test;
  
  test = (int*)malloc(P->wall_no * sizeof(int));
  for(i=0; i<P->wall_no; i++)
    test[i] = 1;
  
  Liste = (int*)malloc(P->wall_no * sizeof(int));
  *Liste_inv = (int*)malloc(P->wall_no * sizeof(int));

  for(i=0; i<P->wall_no; i++){
    if(test[i] != 0){
      Check_mat(P->wall[i]->mat);
      M_inv = mat_inv(P->wall[i]->mat);
      a = 0;
      for(j=i; j<P->wall_no; j++){
	if(a == 0){
	  if(test[j]!=0){
	    s = mat_comp(P->wall[j]->mat,M_inv);

	    if(s == 0){
	      test[i]=0; test[j]=0;
	      Liste[j] = (-1)*(i+1);
	      Liste[i] = (i+1);
	      Liste_inv[0][i] = j;
	      Liste_inv[0][j] = i;
	      a = 1;
	    }
	  }
	}/* if(a==0) */
      }/*for(j)*/
      free_mat(M_inv); M_inv = NULL;
    }
  }/*for(i)*/
  free(test);

  return(Liste);
}/** generate_Liste_with_check(P) **/

/*     s = search_mat(word[i],Liste,P->wall_no); */
static int search_mat(int test, int *Liste, int anz)
{
  int 		i;
  
  i=0; 
  while(i<anz && ((test) != Liste[i]))
    i++;
  /* i ist Position an der test in Liste steht */
  
  return(i);
}/** search_mat() **/

static int search_corner(int s1, int s2, corner_TYP *corners)
{
  int smax, smin;
  int i;
  
  smax = (s1>s2)?s1:s2 ;
  smin = (s1<s2)?s1:s2 ;
  
  i = 0;
  while(smax > corners[i].w[0])
    i++;

  while(smin > corners[i].w[1])
    i++;
  return(i);
}/** search_corner() **/

static int satisfies_map_extra(P,w,vec,M)
     polyeder_TYP *P;
     int w; 
     matrix_TYP *vec; 
     matrix_TYP *M;
{ 
  int i,j,k;
  int e; 
  int dim;
  int *tmp;
  
  dim = P->vert[0]->dim;
  for(i=0; i<P->wall[w]->ext_no; i++){
    e = 0;
    tmp = (int*)calloc(dim,sizeof(int));
    for(j=0; j<dim; j++){
      for(k=0; k<dim; k++)
        tmp[j] += P->wall[w]->extra[i][k] * M->array.SZ[k][j];
    }
    for(j=0; j<dim; j++)
      e += tmp[j] * vec->array.SZ[j][0];
    if(e < 0){
      free(tmp);
      return(-1);
    }
    free(tmp);
  }
  for(i=0; i<P->wall[w]->next_no; i++){
    e = 0;
    tmp = (int*)calloc(dim,sizeof(int));
    for(j=0; j<dim; j++){
      for(k=0; k<dim; k++)
        tmp[j] += P->wall[w]->next[i][k] * M->array.SZ[k][j];
    }
    for(j=0; j<dim; j++)
      e += tmp[j] * vec->array.SZ[j][0];
    if(e < 0){
      free(tmp);
      return(-1);
    }
    free(tmp);
  }
  return(1);
}/** satisfies_map_extra() **/

static int is_on_wall(corner_TYP corn,int *wallungl,polyeder_TYP *P,int d,int w,matrix_TYP *M)
{
  int 	i, j, s, k;
  int 	test; 
  matrix_TYP *v;
  
  v = init_mat(P->vert[0]->dim,1,"");
  test = 0;
  for(j=0; j<corn.v_anz && (test == 0); j++){ 
    for(i=0; i<d; i++)
      test += wallungl[i] * P->vert[corn.v[j]]->v[i];
  }
  if(test == 0){
    for(j=0; j<corn.v_anz && (test == 0); j++){ 
      s = 0;
      for(k=0; k<P->vert[0]->dim; k++)
	v->array.SZ[k][0] = P->vert[corn.v[j]]->v[k];
      s = satisfies_map_extra(P,w,v,M);
      if(s != 1)
	test = 1;
    }
  }
  if(test == 0)
    return(1);
  else
    return(0); 
}/** is_on_wall() **/

static int *map(polyeder_TYP *P, int i, matrix_TYP *M)
{
  int a,b;
  int *erg;
  
  erg = (int*)calloc(M->cols , sizeof(int));
  for(b=0; b< M->cols-1; b++){
    for(a=0; a<M->rows-1; a++){
      erg[b] += P->wall[i]->gl[a] * M->array.SZ[a][b];
    }
  }
  for(a=0; a<M->rows-1; a++){
    erg[M->cols-1] += P->wall[i]->gl[a] * M->array.SZ[a][M->cols-1];
  }
  erg[M->cols-1] += (M->kgv)*(P->wall[i]->gl[M->cols-1]);
  return(erg);
}/* map() */

static int g_wall(int s2, corner_TYP cork,int * Liste_inv, polyeder_TYP *P)
{
  int	i,j;
  int	k;
  int	*test;
  matrix_TYP *M;
  
  /* 3) finde s3 mit k in g2(s3)		*/
  M = P->wall[Liste_inv[s2]]->mat;
  j = 0;
  for(i=0; (i<P->wall_no) && (j == 0); i++){
    if(i != Liste_inv[s2]){
      
      test = map(P,i,M);
      
      j = is_on_wall(cork,test, P, P->vert[0]->dim,i,M);
      /*  cork = corners[k] */
      if( j != 0)
	return(i); 	/* i ist Nr der neuen Wand */
      
      free(test);
    }/* if(i != Liste_inv[s2])*/
  }/* for(i)*/
  printf("Fehler beim Umlaufen der Kante %d \n", k);
  exit(3); 
}/** g_wall() **/

static int create_new_wall(int s2, int s1, int *Liste_inv, int d, polyeder_TYP *P, corner_TYP *corners, int anz)
{ 
  int	s,k;
  int	s3;
  
  /* 1) finde *s1				*/
  s = Liste_inv[s1]; 
  
  /* 2) finde Kante k von *s1 mit s2	*/
  if(s != s2)
    k = search_corner(s,s2,corners);
  else{
    printf("Fehler in create_new_wall ");
    exit(3);
  }
  
  /* 3) finde s3 mit k in g2(s3)		*/
  s3 = g_wall(s2, corners[k],Liste_inv, P);
  
  return(s3);
}/** create_new_wall() **/


/*ACHTUNG REIHENFOLGE -> operiere von links, muss also die Reihenfolge 
  NICHT umdrehen wenn man das "wort" auf einen Punkt anwendet*/
static void append_relation(int S,int *lhsproduct,int lhs_no,int *Liste)
{
  int	ss;
  int 	SIZE;
  
  SIZE = 256;
  if(lhs_no == SIZE-1){
    SIZE += 256;
    lhsproduct = (int *)realloc(lhsproduct, SIZE * sizeof(int)); 
  }
  
  ss = Liste[S-1]; 
  lhsproduct[lhs_no] = ss;
  
}/** append_relation() **/

       
/*		recall(erg->relators[i].lhsproduct[j],erzlist, erg->gen_no);*/
static int recall(int rel, int* erzlist, int gen_no)
{
  int 	j;
  
  if(rel>0){
    for(j=0; j < gen_no; j++){
      if(erzlist[j] == rel)
	return(j+1);
    }
  }
  else{
    for(j=0; j < gen_no; j++){
      if(erzlist[j] == (-1)*rel)
	return(-j-1);
    }
  }
  printf("Fehler in recall \n");
  return(-443);
}/** recall() **/

static matrix_TYP *pres_word_to_mat(int *word, int length, 
				    polyeder_TYP *P, int *Liste)
{
  matrix_TYP 	*A, *B;
  int		i,s;
  
  A = init_mat(P->wall[0]->mat->rows, P->wall[0]->mat->cols,"1");
  
  for(i=0; i<length; i++){
    s = search_mat(word[i],Liste,P->wall_no);
    if(s == P->wall_no){ printf("error in pres_word_to_mat \n"); exit(3);}
    B = mat_mul(A,P->wall[s]->mat);
    Check_mat(B);
    free_mat(A); A = NULL;
    A = B;
  }
  return(A);
}/** pres_word_to_mat() **/

static int calc_order(matrix_TYP *A, matrix_TYP *B)
{ 
  int 	erg;
  int 	count;
  matrix_TYP *test;
  matrix_TYP *tmp;
  
  count = 1;

  test = copy_mat(A);
  while( mat_comp(test,B) != 0){
    tmp = mat_mul(A,test);
    Check_mat(tmp);
    count ++;
    free_mat(test); test = NULL;
    test = copy_mat(tmp);
    free_mat(tmp); tmp = NULL; /*anne*/
  }
  
  erg = count;
  return(erg);
}/** calc_order() **/ 

/**************************************************************************\
@--------------------------------------------------------------------------
@ anne_presentation_TYP *pres()
@
@ polyeder_TYP *Pol      a polyeder with sidetransformations
@
@ calculates a presentation (via a Poincare-algorithm) for the group generated
@ by the sidetransformations of Pol 
@
@--------------------------------------------------------------------------
\**************************************************************************/
anne_presentation_TYP *pres(Pol)
polyeder_TYP *Pol;
{
 int			 anz;
 anne_presentation_TYP 	*erg;
 anne_relator_TYP	*siderel;
 corner_TYP		*corners;
 int 			*Liste, *Liste_inv;
 int 			i,j,s,l,k;
 int 			**Bar;
 matrix_TYP		*bar;
 int 			d;
 int			s1, s2, S1;
 int			neu;
 int			ord;
 int			*erzlist;
 matrix_TYP		*relmat;
 matrix_TYP		*I;
 int			count;
 
 /*Initialisiere*/
 d = Pol->vert[0]->dim;
 bar = init_mat(d,1,"");
 I = init_mat(d,d,"1");
 
 /* erzeuge Liste der Kanten des Fundamental-Polyeders */
 corners = generate_corners(Pol,&anz);
 
 /* Praesentation initialisieren */
 erg = (anne_presentation_TYP *)malloc(sizeof(anne_presentation_TYP));
 erg->generators = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
 erg->relators = (anne_relator_TYP *)malloc(anz * sizeof(anne_relator_TYP));
 erg->rel_no = anz;
 erzlist = (int*)calloc(Pol->wall_no, sizeof(int));
 
 /* Erzeuger abspeichern, dabei Inverse nicht beruecksichtigen*/ 
 
 Liste =  generate_Liste_with_check(Pol,&Liste_inv);
 
 /* Liste der Erzeuger (besteht aus allen Pol->wall[i]->mat)*/
 
 siderel= (anne_relator_TYP *)malloc(anz * sizeof(anne_relator_TYP));
 for(i=0; i<anz; i++){
   siderel[i].lhsproduct = (int*)calloc(2,sizeof(int));
   siderel[i].rhsproduct = (int*)calloc(1,sizeof(int));
   siderel[i].rhsproduct[0] = 0; 
   siderel[i].lhs_no = 2;
   siderel[i].rhs_no = 0;
 }
 count = 0;
 s = 0;
 for(i=0; i<Pol->wall_no; i++){
   if(Liste[i] > 0){
     erg->generators[s] = copy_mat(Pol->wall[i]->mat);
     erzlist[s] = Liste[i];
     if(erzlist[s] != i+1) printf(" erzlist merkwuerdig \n");  /*anne, 30.7.*/
     if(Liste_inv[i] == i){   	/* Seite zu sich selbst gepaart->Relation*/
       siderel[(count)].lhsproduct[0] = s+1; 
       siderel[(count)].lhsproduct[1] = s+1; 
       (count)++;
     }/* if */
     s++;
   }
 }
 erg->gen_no = s;
 erzlist = (int*)realloc(erzlist, s*sizeof(int));
 
siderel = (anne_relator_TYP*)realloc(siderel,(count)*sizeof(anne_relator_TYP));
 
  /*bestimme Baryzentren der Waende von Pol*/
   Bar = bary_list(Pol);

 /* Liste der bereits durchlaufenen Kanten (oder Schleife bis anz=Corner_anz)*/
 /*bestimme Relationen bis auf Ordnung vgl Ratcliffe S.258 */
   for(i=0; i<anz; i++){

   /* fuer jede Kante bestimme Sequenz von Seiten wie folgt:
      D = (w[0],w[1]) = (S1, g1S2) 
      => S1=w[0], 
         S2 erfuellt g1*S2 geschnitten S1 = Corner K
	 S3 erfuellt g2*S3 geschnitten S2 = g1^(-1)*S1 geschn. S2 etc.
      dazu:
      nehme Kante, wird durch waende w[0], w[1] erzeugt. starte mit S1=w[0],
      bilde waende mit (Pol->wall[w[0]]->mat) ab und suche 
      die Wand auf der die Kante liegt => S2.
      Abbilden der Waende : a*x+b=0 => 
      a`*(gx+t)+b` = 0 <=> (a',b') = (a,b)*(g,t)^(-1) = (a*g^(-1),-a*g^(-1)t+b)
   */    

     S1 = corners[i].w[0];
     erg->relators[i].lhsproduct = (int*)calloc(256,sizeof(int));
     erg->relators[i].rhsproduct = (int*)calloc(1,sizeof(int));
     erg->relators[i].rhsproduct[0] = 0;
     erg->relators[i].lhsproduct[0] = Liste[S1]; 
     erg->relators[i].lhs_no = 1;
     erg->relators[i].rhs_no = 0;
     s1 = S1;
     s2 = g_wall(s1, corners[i],Liste_inv, Pol);
     
     l = search_corner(Liste_inv[s1],s2,corners);
     neu = s2;

     while( neu != S1 || i != l)
       {
	 /*anne 18.8. append "neu+1" da Liste von 1 an zaehlt*/
	 append_relation(neu+1,erg->relators[i].lhsproduct,
		         erg->relators[i].lhs_no, Liste);
	 erg->relators[i].lhs_no ++; 
	 
         neu = create_new_wall(s2,s1, Liste_inv, d, Pol, corners, anz);
         s1 = s2;
         s2 = neu; 
	 /*  ist Kante k von *s1 mit s2 =   umlaufende Kante?..*/
	 l = search_corner(Liste_inv[s1],s2,corners);
       }
     /* bestimme noch die Ordnungen der Relatoren */
     for(j=0; j<d; j++)
       bar->array.SZ[j][0] = Bar[corners[i].w[0]][j];   
     l = erg->relators[i].lhs_no;
     relmat = pres_word_to_mat(erg->relators[i].lhsproduct,l,Pol,Liste);
     
     ord = calc_order(relmat,I);
     free_mat(relmat); relmat = NULL; /*anne*/
     
     /*schreibe die Relationen noch in Nr der Erzeuger 
       (wichtig, falls man nur einen Teil der Seitentrafos als Erzeuger hat)*/
     for(j=0; j<l; j++){
       erg->relators[i].lhsproduct[j] = 
	 recall(erg->relators[i].lhsproduct[j],erzlist, erg->gen_no);
     }
     
     /* verlaengere Relation */
     erg->relators[i].lhsproduct = 
       (int *)realloc(erg->relators[i].lhsproduct,
		      (erg->relators[i].lhs_no * ord) * sizeof(int));
     for(j=1; j<ord; j++){
       for(k=0; k<l; k++)
	 erg->relators[i].lhsproduct[j*l +k] = erg->relators[i].lhsproduct[k];
     }
     erg->relators[i].lhs_no = l * ord;
     
   } /* die i-te Kante ist abgearbeitet */    
   free_mat(bar); bar = NULL;
   
   /* nun noch die Seitenpaarungs-relationen hinzufuegen, -> Liste */
   count += erg->rel_no;
   erg->relators = 
    (anne_relator_TYP*)realloc(erg->relators,(count)*sizeof(anne_relator_TYP));
   for(i=erg->rel_no; i<count; i++){
     erg->relators[i].lhsproduct = (int*)calloc(2,sizeof(int));
     erg->relators[i].rhsproduct = (int*)calloc(1,sizeof(int));
     memcpy(&(erg->relators[i]), &(siderel[count-1-i]), sizeof(anne_relator_TYP));
   }
   erg->rel_no = count;
   
   free(erzlist);
   free(Liste);
   free(Liste_inv);
   for(i=0; i<anz; i++)
     free(corners[i].v);
   free(corners);
   
   return(erg);
}/** pres() **/
