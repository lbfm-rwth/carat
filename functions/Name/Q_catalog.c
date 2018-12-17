/****************************************************************************
@
@ ----------------------------------------------------------------------------
@
@ FILE: Q_catalog.c
@
@ ----------------------------------------------------------------------------
@
******************************************************************************/

#include "typedef.h"
#include "name.h"

void (*(display_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data);
void (*(load_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (const char *string, entry *data);
void (*(delete_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data);
int (*(compare_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data1, entry *data2);
const char *name_element [NR_OF_ELEMENTS_IN_EACH_ENTRY];


static int bitfield [32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536}; 

static void *xmalloc(int size, const char *string)
{
  void *pointer;
  if ( (pointer = malloc (size)) == NULL)
    {
      perror (string);
      exit (2);
    }

  return pointer;
}



static void display_abbreviation (entry *element)
{
  fprintf (stdout, "%s  ", element->abbreviation);
}

static void display_degree (entry *element)
{
  fprintf (stdout, "%i  ", element->degree);
}

static void display_symbol (entry *element)
{
  fprintf (stdout, "%s   ", element->symbol);
}


static void display_order (entry *element)
{
  fprintf (stdout, "%i  ", element->order);
}

static void display_discriminant (entry *element)
{
  fprintf (stdout, "%s  ", element->discriminant);
}

static void display_zclasses (entry *element)
{
  fprintf (stdout, "%i  ", element->zclasses);
}

static void display_affine (entry *element)
{
  fprintf (stdout, "%i  ", element->affine);
}

static void display_torsionfree (entry *element)
{
  if (element->torsionfree != -1)
    fprintf (stdout, "%i  ", element->torsionfree);
  else
    fprintf (stdout, "Unknown  ");
}

static void display_no_conclass (entry *element)
{
  fprintf (stdout, "%i   ", element->no_conclass);
}

static void display_no_idem (entry *element)
{
  fprintf (stdout, "%i   ", element->no_idem);
}


static void load_abbreviation (const char *string, entry *element)
{
  /*  if ( (element->abbreviation = (char *) malloc ( (strlen (string) + 1) * sizeof (char))) == NULL     ||
       sprintf (element->abbreviation, string) < 0)
    {
      perror ("load_abbreviation");
      exit (EXIT_FAILURE);
      }*/
  
  element->abbreviation = (char *) xmalloc ( (strlen(string)+1) * sizeof(char),
					     "load_abbreviation");

  if ( sprintf (element->abbreviation, "%s", string) < 0)
    {
      perror ("load_abbreviation");
      exit (4);
    }
  
}

static void load_degree (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->degree)) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_symbol (const char *string, entry *element)
{
  
  /*  if ( (element->symbol = (char *) malloc ( (strlen (string) + 1) * sizeof (char))) == NULL      ||
       sprintf (element->symbol, string) < 0)
    {
      perror ("load_symbol");
      exit (EXIT_FAILURE);
      }*/

  element->symbol = (char *) xmalloc ( (strlen(string)+1) * sizeof(char),
				       "load_symbol");

  if( sprintf(element->symbol, "%s", string) < 0)
    {
      perror ("load_symbol");
      exit (4);
    }
}


static void load_order (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->order)) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_discriminant (const char *string, entry *element)
{
  
/*   if ( (element->discriminant = (char *) malloc ( (strlen (string) + 1) * sizeof (char))) == NULL       || */
/*        sprintf (element->discriminant, string) < 0) */
/*     { */
/*       perror ("load_discriminant"); */
/*       exit (EXIT_FAILURE); */
/*     } */

  element->discriminant = (char *) xmalloc ( (strlen(string)+1) * sizeof(char),
					     "load_discriminant");
  
  if( sprintf (element->discriminant, "%s", string) < 0)
    {
      perror ("load_discriminant");
      exit (4);
    }
  
}

static void load_zclasses (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->zclasses)) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_affine (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->affine)) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_torsionfree (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->torsionfree)) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_no_conclass (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->no_conclass) ) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}

static void load_no_idem (const char *string, entry *element)
{
  if (sscanf (string, "%i", &(element->no_idem) ) != 1)
    {
      fprintf (stderr, "Data has wrong structure.\n");
      exit (4);
    }
}



static void delete_abbreviation (entry *element)
{
  free (element->abbreviation);
}

static void delete_degree (entry *element)
{
  return;
}

static void delete_symbol (entry *element)
{
  free (element->symbol);
}

static void delete_order (entry *element)
{
  return;
}

static void delete_discriminant (entry *element)
{
  free (element->discriminant);
}

static void delete_zclasses (entry *element)
{
  return;
}

static void delete_affine (entry *element)
{
  return;
}

static void delete_torsionfree (entry *element)
{
  return;
}

static void delete_no_conclass (entry *element)
{
  return;
}

static void delete_no_idem (entry *element)
{
  return;
}


static int cmp_abbreviation (entry *entry1, entry *entry2)
{
  /* If the abb(..)-name contains no dot, then it should match with any
     database entry matching to the first dot. */
  if (strchr (entry2->abbreviation, '.') != NULL)
    {
      if (strcmp (entry1->abbreviation, entry2->abbreviation) == 0)
	return 0;
    }
  else
    if (strncmp (entry1->abbreviation, entry2->abbreviation, strlen (entry2->abbreviation)) == 0)
      return 0;
  
  return bitfield [COND_ABBREVIATION];
}

static int cmp_degree (entry *entry1, entry *entry2)
{
  if (entry1->degree == entry2->degree)
    return 0;
  
  return bitfield [COND_DEGREE];
}

static int cmp_symbol (entry *entry1, entry *entry2)
{
  if ( ! strcmp (entry1->symbol, entry2->symbol))
    return 0;
  
  return bitfield [COND_SYMBOL];
}

static int cmp_order (entry *entry1, entry *entry2)
{
  if (entry1->order == entry2->order)
    return 0;
  
  return bitfield [COND_ORDER];
}

static int cmp_discriminant (entry *entry1, entry *entry2)
{
  if ( ! strcmp (entry1->discriminant, entry2->discriminant))
    return 0;
  
  return bitfield [COND_DISCRIMINANT];
}

static int cmp_zclasses (entry *entry1, entry *entry2)
{
  if (entry1->zclasses == entry2->zclasses)
    return 0;
  
  return bitfield [COND_ZCLASSES];
}

static int cmp_affine (entry *entry1, entry *entry2)
{
  if (entry1->affine == entry2->affine)
    return 0;
  
  return bitfield [COND_AFFINE];
}

static int cmp_torsionfree (entry *entry1, entry *entry2)
{
  if (entry1->torsionfree == entry2->torsionfree)
    return 0;
  
  return bitfield [COND_TORSIONFREE];
}

static int cmp_no_conclass (entry *entry1, entry *entry2)
{
  if (entry1->no_conclass == entry2->no_conclass)
    return 0;
  
  return bitfield [COND_NO_CONCLASS];
}

static int cmp_no_idem (entry *entry1, entry *entry2)
{
  if (entry1->no_idem == entry2->no_idem)
    return 0;
  
  return bitfield [COND_NO_IDEM];
}



static void read_database_entry (FILE *file, entry *data)
{
  int i;
  char string[200];
  
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    {
      if (fscanf (file, "%199s", string) != 1)
	{
	  fprintf (stderr, "Data has wrong structure.\n");
	  exit (4);
	}
      (load_element [i]) (string, data);
    }
}


void apply_cond_to_display_list (conditions *cond, database *database, int display_list[], int new_condition)
{
  int i;
  
  for (i=0; i<database->nr; i++)
    {
      display_list [i] = display_list [i] | bitfield [new_condition];
      display_list [i] -= (compare_element [new_condition]) ( & ( (database->entry)[i]), & (cond->entry));
    }
  return;
}

void unapply_cond_to_display_list (database *database, int display_list[], int unset_condition)
{
  int i;
  for (i=0; i<database->nr; i++)
    display_list [i] = display_list [i] | bitfield [unset_condition];
}

database *load_database (const char *filename, int degree)
{
  FILE *file;
  
  int i, j, entries_in_this_file;

  char *complete_name;

  database *datas;

  (display_element [COND_ABBREVIATION]) = display_abbreviation,
    (display_element [COND_DEGREE]) = display_degree,
    (display_element [COND_SYMBOL]) = display_symbol,
    (display_element [COND_ORDER]) = display_order,
    (display_element [COND_DISCRIMINANT]) = display_discriminant,
    (display_element [COND_ZCLASSES]) = display_zclasses,
    (display_element [COND_AFFINE]) = display_affine,
    (display_element [COND_TORSIONFREE]) = display_torsionfree,
    (display_element [COND_NO_CONCLASS]) = display_no_conclass,
    (display_element [COND_NO_IDEM]) = display_no_idem;

  (load_element [COND_ABBREVIATION]) = load_abbreviation,
    (load_element [COND_DEGREE]) = load_degree,
    (load_element [COND_SYMBOL]) = load_symbol,
    (load_element [COND_ORDER]) = load_order,
    (load_element [COND_DISCRIMINANT]) = load_discriminant,
    (load_element [COND_ZCLASSES]) = load_zclasses,
    (load_element [COND_AFFINE]) = load_affine,
    (load_element [COND_TORSIONFREE]) = load_torsionfree,
    (load_element [COND_NO_CONCLASS]) = load_no_conclass,
    (load_element [COND_NO_IDEM]) = load_no_idem;

  (delete_element [COND_ABBREVIATION]) = delete_abbreviation,
    (delete_element [COND_DEGREE]) = delete_degree,
    (delete_element [COND_SYMBOL]) = delete_symbol,
    (delete_element [COND_ORDER]) = delete_order,
    (delete_element [COND_DISCRIMINANT]) = delete_discriminant,
    (delete_element [COND_ZCLASSES]) = delete_zclasses,
    (delete_element [COND_AFFINE]) = delete_affine,
    (delete_element [COND_TORSIONFREE]) = delete_torsionfree,
    (delete_element [COND_NO_CONCLASS]) = delete_no_conclass,
    (delete_element [COND_NO_IDEM]) = delete_no_idem;
  
  (compare_element [COND_ABBREVIATION]) = cmp_abbreviation,
    (compare_element [COND_DEGREE]) = cmp_degree,
    (compare_element [COND_SYMBOL]) = cmp_symbol,
    (compare_element [COND_ORDER]) = cmp_order,
    (compare_element [COND_DISCRIMINANT]) = cmp_discriminant,
    (compare_element [COND_ZCLASSES]) = cmp_zclasses,
    (compare_element [COND_AFFINE]) = cmp_affine,
    (compare_element [COND_TORSIONFREE]) = cmp_torsionfree,
    (compare_element [COND_NO_CONCLASS]) = cmp_no_conclass,
    (compare_element [COND_NO_IDEM]) = cmp_no_idem;
  
  name_element [COND_ABBREVIATION] = "abbreviation",
    name_element [COND_DEGREE] = "degree",
    name_element [COND_SYMBOL] = "symbol",
    name_element [COND_ORDER] = "order",
    name_element [COND_DISCRIMINANT] = "discriminant",
    name_element [COND_ZCLASSES] = "zclasses",
    name_element [COND_AFFINE] = "affine_classes",
    name_element [COND_TORSIONFREE] = "torsionfree_space_groups",
    name_element [COND_NO_CONCLASS] = "con_classes_number",
    name_element [COND_NO_IDEM] = "idempotent_number";


  complete_name = xmalloc ( (strlen(filename) + 16) * sizeof(char), "load_database");
  
  datas = (database *) xmalloc (sizeof (database), "load_database");
  
  datas->nr = 0;

  datas->entry = NULL;

  for (i=0; i < 6; i++)
    if (degree == 0 || degree == i+1)
      {
	
	(void) sprintf (complete_name, "%s%d",filename,i+1);
	
	if ( (file = fopen (complete_name, "r") ) == NULL)
	  {
	    perror ("Open database-file");
	    fprintf (stderr, "Could not open: %s\n", complete_name);
	    exit (4);
	  }
	
	if (fscanf (file, "%i\n", &entries_in_this_file) != 1)
	  {
	    fprintf (stderr, "Data has wrong structure.\n");
	    exit (4);
	  }
	
	if (datas->entry == NULL)
	  /* Dieses if wird durch ein BUG in unseren Malloc-Wrappern
	     notwendig.  Denn "realloc" mit uebergegenem Null-Pointer
	     verhaelt sich so, wie "malloc". Bei unserer Library wird
	     aber, wird 4 vom Nullpointer abgezogen und dann wird das
	     eigentlich "realloc" mit Speicherzelle -4 aufgerufen, was
	     zu einem Segmetation-fault fuehrt. */
	  datas->entry = (entry *) xmalloc (entries_in_this_file * sizeof (entry), "load_database");
	else if ( (datas->entry = (entry *) realloc (datas->entry, (datas->nr + entries_in_this_file) * sizeof (entry)) ) == NULL)
	  {
	    perror ("load database");
	    exit (2);
	  }
	
	for (j=0; j < entries_in_this_file; j++)
	  read_database_entry (file, &(datas->entry [datas->nr + j]) );
	
	datas->nr += entries_in_this_file;
	
	if (fclose (file) == EOF)
	  {
	    perror ("Closing database-file");
	    exit (4);
	  }
	
      }
  
  
  free (complete_name);
  
 



  
  return datas;
}

void free_database (database *datas)
{
  int i, j;

  for (i=0; i<datas->nr; i++)
    for (j=0; j<NR_OF_ELEMENTS_IN_EACH_ENTRY; j++)
	(delete_element [j]) ( &(datas->entry [i]) );


  free (datas->entry);
  free (datas);
}

void display_entry (entry *element)
{
  int i;
  
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    (display_element [i]) (element);
  
  fprintf (stdout, "\n");
}

void display_data_list (database *datas, int display_list[])
{
  int i;
  
  fprintf (stdout,"\n");
  
  for (i=0; i<NR_OF_ELEMENTS_IN_EACH_ENTRY; i++)
    fprintf (stdout,"%s  ", name_element [i]);
  fprintf (stdout,"\n\n");
    

  for (i=0; i<datas->nr; i++)
    if (display_list[i] == ALL_MATCH)
      display_entry ( &(datas->entry[i]) );
  
  fprintf (stdout, "\n");
  
  return;
}
