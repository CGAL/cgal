#include <stdio.h>
#include <string.h>
#include <CGAL/basic.h>
#include "Parse_weight.h"


void usage(char* program)
{
  cerr << "\nNAME\n     "
       << program << " - 2D Alpha shape of a point set\n\n";
  cerr << "SYNOPSIS\n     "
       << program << " [-Delaunay] [-contour] [-regularized]"
       << " [min #] [max #] [winx #] [winy #]"
       << " [-in fname]"
       << " [-out fname]\n";
  cerr << "\nDESCRIPTION\n"
       << "     Computes the alpha shape of a point set that comes from a file or stdin.\n";
  cerr << "\nOPTIONS\n"
       << "     -Delaunay      : Display the underlying Delaunay triangulation" << endl
       << "     -contour       : Display contours" << endl
       << "     -regularized   : Display the regularized alpha shape" << endl
       << "     -min           : xmin and ymin of the logical window" << endl
       << "     -max           : xmax and ymax of the logical window" << endl
       << "     -winx -winy    : size of the physical window" << endl
       << "     -in fname      : Reads points from file ./fname" << endl << endl
       << "     -out fname     : Writes points to file ./fname" << endl << endl
       << "     All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
    strcpy(opt.program, argv[0]);
    --argc;
    argv++;

    while ((argc > 0) && (argv[0][0] == '-')){
        if ((!strcmp(argv[0], "-D")) || (!strcmp(argv[0], "-Delaunay"))) {
            opt.Delaunay = true;
            argv++;
            argc--;
        }
	else if ((!strcmp(argv[0], "-c")) || (!strcmp(argv[0], "-contour"))) {
            opt.contour = true;
            argv++;
            argc--;
        }
        else if ((!strcmp(argv[0], "-r")) || (!strcmp(argv[0], "-regularized"))) {
            opt.regularized = true;
            argv++;
            argc--;
        }
	else if (! strcmp(argv[0], "-min")) {
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

        }else if ((!strcmp(argv[0], "-i")) || (!strcmp(argv[0], "-in"))) {
          strcpy(opt.finname, argv[1]);
          opt.file_input = true;
          argv += 2;
          argc -= 2;
      }else if ((!strcmp(argv[0], "-o")) || (!strcmp(argv[0], "-out"))) {
          strcpy(opt.foutname, argv[1]);
          opt.file_output = true;
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


