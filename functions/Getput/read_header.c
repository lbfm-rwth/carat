#include "typedef.h"
#include "datei.h" // for setup_carat_location
#include "getput.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: read_header.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void read_header(argc,argv)
@ int argc, char *argv[], 
@ 
@   reads the filenames and options for a programm
@   the filenames are stored in:             extern char **FILENAMES;
@   the numer of files is stored in:         extern int *FILEANZ;
@   the options are stored in:               extern char *OPTIONS;
@   the number of options is stored in:      extern int *OPTIONANZ;
@   additional integers to the options in :  extern int *OPTIONNUMBERS;
@   These extern variables are defined in globals.h
@
@   The option are allowed to be alphbetic letters and the function
@   distinguishes between upper and lower letters.
@   In the calling of the program, options have to be set behind a '-'.
@   If one wants read an additional integer 'i' to an option 'p', one has to
@   call this in the form: -p=i
@
@
@    Example:
@         program file1 file2 -p=10 -P -d=-1
@    In this example 
@        FILENAMES[1] = file1
@        FILENAMES[2] = file2
@        FILEANZ = 2
@        OPTIONS[0] = p
@        OPTIONS[1] = P
@        OPTIONS[2] = d
@        OPTIONNUMBERS[0] = 10
@        OPTIONNUMBERS[1] = 0
@        OPTIONNUMBERS[2] = -1
@        OPTIONANZ = 3
@
@     WARNING: a call -CF reads only C as an option, not F.
@              there are no blanks allowed in a word "-p=10"
@     
@
@-------------------------------------------------------------------------
@ int is_option(c)
@ char c;
@
@   The return of this function is 1 if the character c is among the options,
@   otherwise 0.
@   A typical call of this function is:  is_option('p');
@
@-------------------------------------------------------------------------
@
@ int optionnumber(c)
@ char c;
@
@   optionnumber('p') returns the additional number to the option 'p'.
@   If 'p' is no option, the return is 0.
@-------------------------------------------------------------------------
@
\**************************************************************************/

char **FILENAMES;
int FILEANZ;
static char *OPTIONS;
static int *OPTIONNUMBERS;
static int OPTIONANZ;

void 
read_header (int argc, char *argv[])
{
  int i;
  char *w;

  FILENAMES = (char **)malloc(argc *sizeof(char *));
  FILEANZ = 0;
  OPTIONS = (char *)malloc(argc *sizeof(char));
  OPTIONNUMBERS = (int *)malloc(argc *sizeof(int));
  OPTIONANZ = 0;
  for ( i = 1; i < argc; i++)
  {
     switch ( argv[i][0] )
     {
        case '-' :
          OPTIONS[OPTIONANZ] = argv[i][1];
          if ( (w = strchr (argv[i], '=')) != NULL )
                    sscanf(w, "=%d", &OPTIONNUMBERS[OPTIONANZ]);
          else
            OPTIONNUMBERS[OPTIONANZ] = 0;
          OPTIONANZ++;
        break;
        default  :
           FILENAMES[FILEANZ] = argv[i];
           FILEANZ++;
     }
  }
  setup_carat_location(argv[0]);
} 


int 
is_option (char c)
{
  int i;
  for(i=0;i<OPTIONANZ;i++)
  {
    if(OPTIONS[i] == c)
      return(TRUE);
  }
  return(FALSE);
}


int 
optionnumber (char c)
{
  int i;
  for(i=0;i<OPTIONANZ;i++)
  {
    if(OPTIONS[i] == c)
      return(OPTIONNUMBERS[i]);
  }
  return(0);
}

void 
unread_header (void)
{
   free(FILENAMES);
   free(OPTIONS);
   free(OPTIONNUMBERS);
}
