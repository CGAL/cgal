#include <CGAL/basic.h>
#include <iostream>
#include <cstring>
#include <cstdio>
#include "parse.h"


void usage(char* program)
{
  std::cerr << "\nNAME\n     "
	    << program << " - Triangulation of a point set\n\n";
  std::cerr << "SYNOPSIS\n     "
	    << program << " [-draw] [-statistics] [-check]"
	    << " [min #] [max #] [winx #] [winy #]"
	    << " [-file fname]\n";

  std::cerr << "\nDESCRIPTION\n"
	    << " Triangulates a point set that comes from a file or stdin.\n";
  std::cerr << "\nOPTIONS\n"
	    << " -draw          : Displays intermediate results" << std::endl
	    << " -statistics    : Collects and displays performance data" 
	    << std::endl
	    << " -check         : Performs correctness tests" << std::endl
	    << " -min           : xmin and ymin of the logical window" 
	    << std::endl
	    << " -max           : xmax and ymax of the logical window" 
  	    << std::endl
	    << "     -winx -winy    : size of the physical window" << std::endl
	    << "     -file fname    : Reads points from file ./fname" 
	    << std::endl << std::endl
	    << "All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
    CGAL_CLIB_STD::strcpy(opt.program, argv[0]);
    --argc;
    argv++;

    while ((argc > 0) && (argv[0][0] == '-')){
        if ((!CGAL_CLIB_STD::strcmp(argv[0], "-d")) || 
	    (!CGAL_CLIB_STD::strcmp(argv[0], "-draw"))) {
            opt.draw = true;
            argv++;
            argc--;
        }
        else if ((!CGAL_CLIB_STD::strcmp(argv[0], "-s") || 
		  (!CGAL_CLIB_STD::strcmp(argv[0], "-statistics")))){
            opt.statistics = true;
            argv++;
            argc--;
        }
        else if ((!CGAL_CLIB_STD::strcmp(argv[0], "-c")) || 
		 (!CGAL_CLIB_STD::strcmp(argv[0], "-check"))) {
            opt.check = true;
            argv++;
            argc--;
        }else if (! CGAL_CLIB_STD::strcmp(argv[0], "-min")) {
            if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.min) != 1) {
                std::cerr << "Argument for min must be a number"
                     << std::endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! CGAL_CLIB_STD::strcmp(argv[0], "-max")) {
            if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.max) != 1) {
                std::cerr << "Argument for max must be a number"
                     << std::endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! CGAL_CLIB_STD::strcmp(argv[0], "-winx")) {
            if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.winx) != 1) {
                std::cerr << "Argument for winx must be a number"
                     << std::endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! CGAL_CLIB_STD::strcmp(argv[0], "-winy")) {
            if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.winy) != 1) {
                std::cerr << "Argument for winy must be a number"
                     << std::endl;
            }
            argv += 2;
            argc -= 2;

        }else if ((!CGAL_CLIB_STD::strcmp(argv[0], "-f")) || 
		  (!CGAL_CLIB_STD::strcmp(argv[0], "-file"))) {
          CGAL_CLIB_STD::strcpy(opt.fname, argv[1]);
          opt.file_input = true;
          argv += 2;
          argc -= 2;
      }
      else if ((!CGAL_CLIB_STD::strcmp(argv[0], "-?")) ||
               (!CGAL_CLIB_STD::strcmp(argv[0], "-h")) ||
               (!CGAL_CLIB_STD::strcmp(argv[0], "-help"))) {
          usage(opt.program);
          return false;
      }
      else {
          std::cerr << "Unrecognized option " << argv[0] << std::endl;
          usage(opt.program);
          return false;
      }
    }
  if(argc > 0){
      std::cerr << "Unrecognized option " << argv[0] << std::endl;
      usage(opt.program);
      return false;
  }
  return true;
}


