/* 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   PRC, the profile comparer, version 1.5.6

	convert_to_prc.c: wrapper for converter to the PRC file format

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2003-5 Martin Madera and MRC LMB, Cambridge, UK
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public 
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prc.h"


static char usage[] =
"Usage: convert_to_prc [-human] <input file> [<output file>]\n"
"\n"
"Recognized input file extensions:\n"
"\n"
KNOWN_EXTENSIONS
"\n"
"The default output format is the current PRC binary format. When\n"
"the '-human' option is used, the output will be a human-readable\n"
"ASCII format similar to SAM. The '-human' format cannot be read\n"
"by PRC, but it may be useful for debugging and experimentation.\n"
"\n"
"When no output file is specified, the output is to STDOUT.\n"
"\n";


// need to have usage[] for common.c
#include "common.c"


int main(int argc, char **argv)
{ 
  FILE  *file;
  int   argi, human = 0;
  HMM   *hmm;


  sprintf( version, "CONVERT-TO-PRC " VERSION " (%s), compiled on %s",
#if   PROF_HMM_TRANS==PLAN9
	   "PLAN9"
#elif PROF_HMM_TRANS==PLAN7
	   "PLAN7"
#endif
	   , __DATE__);

  // print out the header
  print_header(stderr, "");

  // sort out the options
  for( argi=1; argi<argc; argi++ )
    {
      // onto the filenames ...
      if( argv[argi][0]!='-' ) break;

      if(strncmp(argv[argi],"-human",7)==0)
	{
	  human = 1;
	}
     else if( (strncmp(argv[argi],"-h",3)==0) ||
	       (strncmp(argv[argi],"--help",7)==0) )
	{
	  print_help_die();
	}
      else
	{
	  arg_error("Error: Unknown option '%s'!", argv[argi]);
	};      
    };

  // check that have one or two filenames
  if( argc-argi == 0 )
    print_help_die();
  else if( argc-argi > 2 )
    arg_error("Error: Incorrect number of arguments!");
  

  hmm = read_HMM(argv[argi++]);

  if( argc-argi == 1 )
    {
      open_file_or_die(file, argv[argi], "w");
    }
  else
    {
      file = stdout;
    };

  if( human )
    print_HMM(file, hmm);
  else
    write_HMM_PRC_binary(file, hmm);

  if( argc-argi == 1 )
    fclose(file);

  free_HMM(hmm);

  return 0;
}
