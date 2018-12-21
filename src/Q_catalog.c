#include "name.h"
#include "getput.h"
#include "sort.h"
#include "bravais.h"
#include "gmp.h"
#include "name.h"
#include "matrix.h"
#include "datei.h"

int SFLAG;
int INFO_LEVEL;

static void write_groups_to_file (int *write_list, database *database, char *string)
{
  char argument [80], inputfilename [1024], format[1024];
  FILE *outputfile, *inputfile;
  int i, data_char;
  
  if (sscanf (string, "w%s", argument) == 1)
    {
      if (strcmp (argument, "-s"))
	{
	  if ( (outputfile = fopen (argument, "w") ) == NULL)
	    {
	      perror ("Open output-file");
	      fprintf (stderr, "Could not open: %s./n", argument);
	      exit (EXIT_FAILURE);
	    }	  
	}
      else
	{
	  printf("stdout\n");
	  outputfile = stdout;
	}
    }
  else
    {
      fprintf (stdout, "Write to file: w <filename>\nWrite to stdout: w -s\n\n\n");
      return;
    }

  get_data_dir(format, "tables/qcatalog/dim%i/dir.%s/ordnung.%i/%s/%s");
  for (i=0; i<database->nr; i++)
    if (write_list[i] == ALL_MATCH)
      {
	sprintf (inputfilename, format, 
		 (database->entry [i]).degree,
		 (database->entry [i]).symbol,
		 (database->entry [i]).order,
		 (database->entry [i]).discriminant,
		 (database->entry [i]).abbreviation);

	
	if ( (inputfile = fopen (inputfilename, "r") ) == NULL)
	  {
	    perror ("Open group-file");
	    fprintf (stderr, "Could not open: %s./n", inputfilename);
	    exit (EXIT_FAILURE);
	  }

	while ( (data_char = fgetc (inputfile)) != EOF)
	  {
	    fputc (data_char, outputfile);
	  }

	if (fclose (inputfile) == EOF)
	  {
	    perror ("Closing input-group-file");
	    exit (EXIT_FAILURE);
	  }
  
  
      }
  
  if (outputfile != stdout)
    if (fclose (outputfile) == EOF)
      {
	perror ("Closing output-file");
	exit (EXIT_FAILURE);
      }
  
  return;
}

void display_conditions (conditions *cond)
{
  int i;
  
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    if ( (cond->exists) [i] == TRUE)
      {
	fprintf (stdout, " %s = ", name_element [i]);
	(display_element [i]) (&(cond->entry));
      }
}

/* This function tests how many possible arguments of type "sub_entry"
   still exist with the choice of other conditions
   given. "list_of_possible" returns the list of the elements still
   possible given by the array number to look it up in the "database"
   and the return value of this function is the number of possible
   values. */
int possible_arguments (int *search_list, database *database,
			int **list_of_possible, int sub_entry)
{
  int i, j, counter = 0;

  char  *new_element;
  
  new_element = (char *) malloc(database->nr * sizeof(char));
  
  new_element = memset (new_element, TRUE, database->nr);


  (*list_of_possible) = (int *) malloc(database->nr * sizeof(int));
  
  for (i=0; i<database->nr; i++)
    if (search_list [i] == ALL_MATCH  &&  new_element[i] == TRUE)
      {
	(*list_of_possible)[counter] = i;
	counter ++;
	
	for (j=i+1; j< database->nr; j++)
	  if (search_list [j] == ALL_MATCH  &&  new_element[j] == TRUE  &&
	      (compare_element [sub_entry])(  &(database->entry[i]),
					      &(database->entry[j]))
	      == 0)
	    new_element [j] = 0;
      }

  free (new_element);

  (*list_of_possible) = (int *) realloc( (*list_of_possible), counter * sizeof(int));

  return counter;
}

void display_info_text (int *search_list, database *database,
			conditions *cond)
{
  int i,
    q_counter = 0,
    z_counter = 0,
    aff_counter = 0;
  
  char part1[4];
  
  fprintf (stdout, "Conditions:  ");
  display_conditions (cond);


  for (i=0; i<database->nr; i++)
    if (search_list [i] == ALL_MATCH){
      q_counter ++; 
      z_counter += database->entry[i].zclasses;
      aff_counter += database->entry[i].affine;
    }
 
  fprintf (stdout, "\nQ_Classes: # %i,   Z_Classes: # %i,   Affine Classes: #%i\n\n",q_counter, z_counter, aff_counter);
  


  
  fprintf (stdout, "Possible Conditions:  ");
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    if (strlen (name_element [i]) > 3)
      {
	part1[0] = (name_element [i]) [0],
	  part1[1] = (name_element [i]) [1],
	  part1[2] = (name_element [i]) [2],
	  part1[3] = '\00';
	

	fprintf (stdout, "%3s(%s)   ", part1, &(name_element [i][3]));
      }
    else
      fprintf (stdout, "%s   ", name_element [i]);
  fprintf (stdout, "\n");

  fprintf (stdout, "\n\ns   set condition		     p     possible data\n");
  fprintf (stdout, "d   delete condition		     l     list group\n");
  fprintf (stdout, "w   write to file                    q     quit\n\n\n");
  
}

int change_conditions (conditions *cond, char *string, char set_or_del)
{
  char part1 [80], part2 [80];
  
  int i;
  
  if (set_or_del == SET_COND)
    {


      if (sscanf (string, "s%s %s", part1, part2) != 2)
	{
	  fprintf (stdout, "Wrong syntax. Please use:\n     s <cond> <value>\n\n\n");
	  return -1;
	}
      
      for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
	if ( (strlen (part1) == 3   &&  strncmp (part1, name_element [i], 3) == 0)
	     || strcmp (part1, name_element [i]) == 0)
	  {
	    if ( (cond->exists) [i] == TRUE )
	      (delete_element [i]) ( &(cond->entry) );
	    (cond->exists) [i] = TRUE;
	    (load_element [i]) (part2, &(cond->entry) );
	    fprintf (stdout, "'%s' has been set to ", name_element [i]);
	    (display_element [i]) ( &(cond->entry) );
	    fprintf (stdout, ".\n\n");
	    return i;
	  }	
      
    }
  else if (set_or_del == DEL_COND)
    {

      sscanf (string, "d%s", part1);

      for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
	if ( (strlen (part1) == 3   &&  strncmp (part1, name_element [i], 3) == 0)
	     || strcmp (part1, name_element [i]) == 0)
	  if ((cond->exists) [i] == TRUE)
	    {
	      (delete_element [i]) ( &(cond->entry) );
	      fprintf (stdout, "'%s' has been unset..\n\n\n", name_element [i]);
	      (cond->exists) [i] = FALSE;
	      return i;
	    }
	  else 
	    return -1;
      
    }
  else if (set_or_del == DISPLAY_POSSIBLE)
    {
      
      sscanf (string, "p%s", part1);
      
      for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
	if ( (strlen (part1) == 3   &&  strncmp (part1, name_element [i], 3) == 0)
	     || strcmp (part1, name_element [i]) == 0)
	  return i;
	    
      
    }
  
  fprintf (stdout, "\n\nWrong input!\n\n");
  return -1;
}

int prompt_input (int *search_list, database *database, conditions *cond)
{
  char string [80], part1 [80];
  int i;

  fprintf (stdout, "\n>");
  for (i=0; i<80; i++)
    if ( (string [i] = fgetc(stdin)) == '\n')
      {
	string [i] = '\000';
	break;
      }

  sscanf (string, "%s", part1);

  switch (part1[0])
    {
    case 'q':
      return 1;
    case 's':
      if ( (i = change_conditions (cond, string, SET_COND)) == -1)
	return 0;
      else
	apply_cond_to_display_list (cond, database, search_list, i);
      break;
    case 'd':
      if ( (i = change_conditions (cond, string, DEL_COND)) == -1)
	return 0;
      else
	unapply_cond_to_display_list (database, search_list, i);
      break;
    case 'p':
      if ( (i = change_conditions (cond, string, DISPLAY_POSSIBLE)) == -1)
	return 0;
      else
	{
	  int *list_of_possible;
	  int number_of_possible, j;

	  number_of_possible = possible_arguments  (search_list, database,
						    &list_of_possible, i);

	  fprintf (stdout, "Possible input for %s:\n", name_element [i]);
	  for (j=0; j<number_of_possible; j++){
	    if ( ! (j % 4) )
	       fprintf (stdout, "\n");
            else
               fprintf (stdout, "\t");
	    (display_element [i]) ( & (database->entry[list_of_possible[j]]) );
	  }

	  fprintf (stdout, "\n\n");


	  free (list_of_possible);
	  
	}
      break;
    case 'l':
      display_data_list (database, search_list);
      return 0;
    case 'w':
      write_groups_to_file (search_list, database, string);
      return 0;
    default:
      fprintf (stdout, "Command '%c' is not defined.\n", part1[0]);
      return 0;
    }
  
  return 0;
}

void interactive_mode (int *search_list, database *database)
{
  conditions *cond;
  int i;
  
  
  if ( (cond = (conditions *) malloc (sizeof (conditions)) ) == NULL  ||
       (cond->exists = (int *) malloc (NR_OF_ELEMENTS_IN_EACH_ENTRY * 
				       sizeof (int)) ) == NULL )
    {
      perror ("interactive_mode");
      exit (EXIT_FAILURE);
    }
  
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    (cond->exists) [i] = FALSE;
  
  
  while (1)
    {
      display_info_text (search_list, database, cond);
      if (prompt_input (search_list, database, cond) != 0)
	{
	  fprintf (stdout, "Program aborted\n");
	  break;
	}
    }
  return;
}

int main (int argc, char *argv[])
{
  int i, *display_list;
  
  database *database;
  
  bravais_TYP *G,
              *H;

  matrix_TYP *T;

  char name[1024], symb[1024], dbname[1024];

  extern int FILEANZ;
  extern char **FILENAMES;

  read_header (argc, argv);
  
  if (is_option('h'))
    INFO_LEVEL = optionnumber('h');
  
  if (INFO_LEVEL == 8)
    SFLAG = 1; 
  
  if (is_option('h') && INFO_LEVEL != 8)
    {
      printf("Usage: %s [file] [-T] [-h] [-i] [-s]\n",argv[0]);
      printf("\n");
      printf("file: bravais_TYP containing the finite unimodular group G.\n");
      printf("\n");
      printf("The program Q_catalog as two identities:\n");
      printf("If called without an input file, it gives access to the\n");
      printf("database of all Q-classes of finite unimodular groups of degree\n");
      printf("up to 6. In this mode, the command 'h' will provide further help.\n");
      printf("If called with an input file, it searches for the given group\n");
      printf("in the database, and gives a unique name choosen for this\n");
      printf("Q-class.\n");
      printf("\n");
      printf("Options:\n");
      printf("-h    : gives this help.\n");
      printf("-T    : calculate a transformation matrix transforming the input\n");
      printf("        group G into the group given in the catalog. (2nd mode only).\n");
      printf("-i    : output the group in the catalog which is Q-equivalent\n");
      printf("        to the input group G. (2nd mode only).\n");
      printf("-s    : output the family symbol of G. (2nd mode only).\n");
      printf("\n");
      printf("Cf: Symbol, Bravais_type, Q_equiv, Bravais_catalog, Conj_bravais.\n");
      exit (EXIT_SUCCESS);
    }
  
  get_data_dir(dbname, "tables/qcatalog/data");
  if (FILEANZ == 1){
    G = get_bravais(FILENAMES[0]);

    database = load_database (dbname, G->dim);
  
    display_list = (int *) malloc (database->nr * sizeof (int));

    for (i=0; i<database->nr; i++)
      display_list [i] = ALL_MATCH;
  
    if (is_option('i')){
      T = q_class_inf (G, database, name, symb, &H, NULL, is_option('T'));
    } 
    else { 
      T = q_class_inf (G, database, name, symb, NULL, NULL, is_option('T'));
    }

    printf("Name of this Q-class: %s\n",name);
    if (is_option('s'))
       printf("symbol of the group %s\n",symb);
    if (is_option('T')){
       put_mat(T,NULL,"transformation matrix",0);
       free_mat(T);
    }
    if (is_option('i')){
       put_bravais(H,NULL,NULL);
       free_bravais(H);
    }

    free_bravais(G);

  }
  else if (FILEANZ == 0){
    database = load_database (dbname, 0);
    
    display_list = (int *) malloc (database->nr * sizeof (int));
    
    for (i=0; i<database->nr; i++)
      display_list [i] = ALL_MATCH;
    
    interactive_mode (display_list, database);
  }
  else
    {
      fprintf (stderr, "Wrong number of files. Use <%s -h> for help.\n", argv[0]);
      exit (EXIT_FAILURE);
    } 
  /*  display_data_list (database, display_list); */

  free(display_list);  
  free_database (database);
  if (INFO_LEVEL == 8) pointer_statistics(0,0);
  
  return 0;
}
