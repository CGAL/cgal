#include <cstdio>
#include <cstring>
#include <CGAL/basic.h>
#include "Parse_weight.h"


void usage(char* program)
{
  std::cerr << "\nNAME\n     "
	    << program << " - 2D Alpha shape of a point set\n\n";
  std::cerr << "SYNOPSIS\n     "
	    << program << " [-Delaunay] [-contour] [-regularized]"
	    << " [min #] [max #] [winx #] [winy #]"
	    << " [-in fname]"
	    << " [-out fname]\n";
  std::cerr << "\nDESCRIPTION\n"
	    << "     Computes the alpha shape of a point set that comes from a file or stdin.\n";
  std::cerr << "\nOPTIONS\n"
	    << "     -Delaunay      : Display the underlying Delaunay triangulation" << std::endl
	    << "     -contour       : Display contours" << std::endl
	    << "     -regularized   : Display the regularized alpha shape" << std::endl
	    << "     -min           : xmin and ymin of the logical window" << std::endl
	    << "     -max           : xmax and ymax of the logical window" << std::endl
	    << "     -winx -winy    : size of the physical window" << std::endl
	    << "     -in fname      : Reads points from file ./fname" << std::endl << std::endl
	    << "     -out fname     : Writes points to file ./fname" << std::endl << std::endl
	    << "     All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
  std::strcpy(opt.program, argv[0]);
  --argc;
  argv++;

  while ((argc > 0) && (argv[0][0] == '-')){
    if ((!std::strcmp(argv[0], "-D")) || (!std::strcmp(argv[0], "-Delaunay"))) {
      opt.Delaunay = true;
      argv++;
      argc--;
    }
    else if ((!std::strcmp(argv[0], "-c")) || (!std::strcmp(argv[0], "-contour"))) {
      opt.contour = true;
      argv++;
      argc--;
    }
    else if ((!std::strcmp(argv[0], "-r")) || (!std::strcmp(argv[0], "-regularized"))) {
      opt.regularized = true;
      argv++;
      argc--;
    }
    else if (! std::strcmp(argv[0], "-min")) {
      if (std::sscanf(argv[1], "%lf", &opt.min) != 1) {
	std::cerr << "Argument for min must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;

    }else if (! std::strcmp(argv[0], "-max")) {
      if (std::sscanf(argv[1], "%lf", &opt.max) != 1) {
	std::cerr << "Argument for max must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;

    }else if (! std::strcmp(argv[0], "-winx")) {
      if (std::sscanf(argv[1], "%d", &opt.winx) != 1) {
	std::cerr << "Argument for winx must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;

    }else if (! std::strcmp(argv[0], "-winy")) {
      if (std::sscanf(argv[1], "%d", &opt.winy) != 1) {
	std::cerr << "Argument for winy must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;

    }else if ((!std::strcmp(argv[0], "-i")) || (!std::strcmp(argv[0], "-in"))) {
      std::strcpy(opt.finname, argv[1]);
      opt.file_input = true;
      argv += 2;
      argc -= 2;
    }else if ((!std::strcmp(argv[0], "-o")) || (!std::strcmp(argv[0], "-out"))) {
      std::strcpy(opt.foutname, argv[1]);
      opt.file_output = true;
      argv += 2;
      argc -= 2;
    }
    else if ((!std::strcmp(argv[0], "-?")) ||
	     (!std::strcmp(argv[0], "-h")) ||
	     (!std::strcmp(argv[0], "-help"))) {
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


