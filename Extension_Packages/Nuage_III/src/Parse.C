#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>

#include <NUAGE/Parse.h>


void usage(char* program)
{
  std::cerr << std::endl << "NAME     " << std::endl
	    << program << "     - surface extension -" << std::endl << std::endl;

  std::cerr << std::endl << "OPTIONS" << std::endl
	    << "     -Delaunay      : display the underlying Delaunay triangulation" << std::endl
	    << "     -contours      : display contours" << std::endl
	    << "     -no_border -nb : set K = infinity" << std::endl
	    << "     -delta x       : set the delta constant" << std::endl
	    << "     -ki x y        : set the K interval (default : [1.1  5])" << std::endl
	    << "     -ks x          : set the K step (default : .1)" << std::endl
	    << "     -k  x          : set the K constant (only one pass)" << std::endl
	    << "     -number_of_points x                : set a number of points for a sub-sample" << std::endl
	    << "     -max_of_connected_components x     : set the max of connected components" << std::endl
	    << "                                          (default : non-active)" << std::endl
            << "     -post x        : set a number for the post process" << std::endl
	    << "     -in fname      : reads points from file ./fname" << std::endl
	    << "     -sect_in fname : reads points from sections file ./fname" << std::endl
	    << "     -out fname     : writes points to file ./fname" << std::endl
	    << "     -out_format -of: choose file format for output (iv, vrml, oogl, medit," << std::endl
	    << "                                                     ply, all, none)" << std::endl
	    << "     All options can be abbreviated by their first character" << std::endl;
}



bool
parse(int argc, char* argv[], Options &opt)
{
  AF_CGAL_CLIB_STD::strcpy(opt.program, argv[0]);
  --argc;
  argv++;
  if(argc == 0)
    std::cout << "nothing ???" << std::endl;
  
  while ((argc > 0) && (argv[0][0] == '-')){
    if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-D")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-Delaunay"))) {
      opt.Delaunay = true;
      argv++;
      argc--;
      std::cout << "-D ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-c")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-contours"))) {
      opt.contour = true;
      argv++;
      argc--;
      std::cout << "-c ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-nb")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-no_border"))) {
      opt.K = HUGE_VAL;
      opt.K_init = opt.K;
      argv++;
      argc--;
      std::cout << "-nb ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-d")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-delta"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.DELTA) != 1) {
	std::cerr << "Argument for DELTA must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-d " << opt.DELTA << " "; 
    }  
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-ki"))){
      if ((CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K_init) != 1)||
	  (CGAL_CLIB_STD::sscanf(argv[2], "%lf", &opt.K) != 1)){
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 3;
      argc -= 3;
      std::cout << "-ki " << opt.K_init << " " << opt.K << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-ks"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K_step) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-ks " << opt.K_step << " ";
    }  
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-k"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      opt.K_init = opt.K;
      argv += 2;
      argc -= 2;
      std::cout << "-k " << opt.K_init << " ";
    } 
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-n")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-number_of_points"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.number_of_points) != 1) {
	std::cerr << "Argument for the number of points must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-n " << opt.number_of_points << " ";
    }  
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-m")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-max_of_connected_components"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.max_connected_comp) != 1) {
	std::cerr << "Argument for the number of connected components must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-m " << opt.max_connected_comp << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-p")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-post"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.NB_BORDER_MAX) != 1) {
	std::cerr << "Argument for post process must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-p " << opt.NB_BORDER_MAX << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-i")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-in"))) {
      AF_CGAL_CLIB_STD::strcpy(opt.finname, argv[1]);
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cout << "-i " << opt.finname << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-s")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-sect_in"))) {
      AF_CGAL_CLIB_STD::strcpy(opt.finname, argv[1]);
      opt.Section_file = true;
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cout << "-s " << opt.finname << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-o")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-out"))) {
      AF_CGAL_CLIB_STD::strcpy(opt.foutname, argv[1]);
      opt.file_output = true;
      argv += 2;
      argc -= 2;
      std::cout << "-o " << opt.foutname << " ";
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-of")) || (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-out_format"))) {
      if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "vrml"))
	opt.out_format = 0;
      else if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "oogl"))
	opt.out_format = 1;
      else if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "medit"))
	opt.out_format = 2;
      else if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "ply"))
	opt.out_format = 3;
      else if(!AF_CGAL_CLIB_STD::strcmp(argv[1], "iv"))
	opt.out_format = 4;
      else if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "all"))
	opt.out_format = -1;
      else if (!AF_CGAL_CLIB_STD::strcmp(argv[1], "none"))
	opt.out_format = -2;
      else
	std::cout << "unrecognized file format." << std::endl;
      opt.file_output = true;
      std::cout << "-of " << argv[1] << " ";
      argv += 2;
      argc -= 2;
    }
    else if ((!AF_CGAL_CLIB_STD::strcmp(argv[0], "-?")) ||
	     (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-h")) ||
	     (!AF_CGAL_CLIB_STD::strcmp(argv[0], "-help"))) {
      usage(opt.program);
      return false;
    }
    else {
      std::cerr << "unrecognized option " << argv[0] << std::endl;
      usage(opt.program);
      return false;
    }
  }
  if(argc > 0){
    std::cerr << "unrecognized option " << argv[0] << std::endl;
    usage(opt.program);
    return false;
  }
  return true;
}


