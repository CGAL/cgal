#include <cstdio>
#include <cstring>
#include <CGAL/basic.h>
#include "parse.h"


void usage(char* program)
{
  cerr << "\nNAME\n     "
       << program << " - Triangulation of a point set\n\n";
  cerr << "SYNOPSIS\n     "
       << program << " [-draw] [-statistics] [-check]"
       << " [min #] [max #] [winx #] [winy #]"
       << " [-file fname]\n";

  cerr << "\nDESCRIPTION\n"
       << "     Triangulates a point set that comes from a file or stdin.\n";
  cerr << "\nOPTIONS\n"
       << "     -draw          : Displays intermediate results" << endl
       << "     -statistics    : Collects and displays performance data" << endl
       << "     -check         : Performs correctness tests" << endl
       << "     -min           : xmin and ymin of the logical window" << endl
       << "     -max           : xmax and ymax of the logical window" << endl
       << "     -winx -winy    : size of the physical window" << endl
       << "     -file fname    : Reads points from file ./fname" << endl << endl
       << "     All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
    strcpy(opt.program, argv[0]);
    --argc;
    argv++;

    while ((argc > 0) && (argv[0][0] == '-')){
        if ((!strcmp(argv[0], "-d")) || (!strcmp(argv[0], "-draw"))) {
            opt.draw = true;
            argv++;
            argc--;
        }
        else if ((!strcmp(argv[0], "-s") || (!strcmp(argv[0], "-statistics")))){
            opt.statistics = true;
            argv++;
            argc--;
        }
        else if ((!strcmp(argv[0], "-c")) || (!strcmp(argv[0], "-check"))) {
            opt.check = true;
            argv++;
            argc--;
        }else if (! strcmp(argv[0], "-min")) {
            if (sscanf(argv[1], "%lf", &opt.min) != 1) {
                cerr << "Argument for min must be a number"
                     << endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! strcmp(argv[0], "-max")) {
            if (sscanf(argv[1], "%lf", &opt.max) != 1) {
                cerr << "Argument for max must be a number"
                     << endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! strcmp(argv[0], "-winx")) {
            if (sscanf(argv[1], "%d", &opt.winx) != 1) {
                cerr << "Argument for winx must be a number"
                     << endl;
            }
            argv += 2;
            argc -= 2;

        }else if (! strcmp(argv[0], "-winy")) {
            if (sscanf(argv[1], "%d", &opt.winy) != 1) {
                cerr << "Argument for winy must be a number"
                     << endl;
            }
            argv += 2;
            argc -= 2;

        }else if ((!strcmp(argv[0], "-f")) || (!strcmp(argv[0], "-file"))) {
          strcpy(opt.fname, argv[1]);
          opt.file_input = true;
          argv += 2;
          argc -= 2;
      }
      else if ((!strcmp(argv[0], "-?")) ||
               (!strcmp(argv[0], "-h")) ||
               (!strcmp(argv[0], "-help"))) {
          usage(opt.program);
          return false;
      }
      else {
          cerr << "Unrecognized option " << argv[0] << endl;
          usage(opt.program);
          return false;
      }
    }
  if(argc > 0){
      cerr << "Unrecognized option " << argv[0] << endl;
      usage(opt.program);
      return false;
  }
  return true;
}


