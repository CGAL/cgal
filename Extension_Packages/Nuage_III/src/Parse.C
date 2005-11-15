#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>

#include <NUAGE/Parse.h>

namespace NUAGE {

void usage(char* program)
{
  std::cerr << std::endl << "NAME     " << std::endl
	    << program << "     - surface extension -" << std::endl << std::endl;

  std::cerr << std::endl << "OPTIONS" << std::endl
	    << "     -Delaunay      : display the underlying Delaunay triangulation" << std::endl
	    << "     -contours      : display contours" << std::endl
	    << "     -shuffle       : random shuffle" << std::endl
	    << "     -binary        : binary I/O" << std::endl
	    << "     -xyz           : input data in xyz format" << std::endl
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
	    << "     -out_format -of: choose file format for output (iv, wrl, oogl, medit," << std::endl
	    << "                                                     ply, stl, all, none)" << std::endl
            << "     -rgb r g b     : color of the surface" << std::endl
            << "     -no_header     : The Vrml header and footer are not written" << std::endl	    
	    << "     -area  a       : set the area threshold: No faces larger than area * average_area" << std::endl
	    << "     -perimeter  p  : set the perimeter threshold: No faces larger than perometer * average_perimeter" << std::endl
	    << "     All options can be abbreviated by their first character" << std::endl;
}



bool
parse(int argc, char* argv[], Options &opt)
{
  std::strcpy(opt.program, argv[0]);
  --argc;
  argv++;
  if(argc == 0)
    std::cout << "nothing ???" << std::endl;
  
  while ((argc > 0) && (argv[0][0] == '-')){
    if ((!std::strcmp(argv[0], "-D")) || (!std::strcmp(argv[0], "-Delaunay"))) {
      opt.Delaunay = true;
      argv++;
      argc--;
      std::cout << "-D ";
    }
    else if ((!std::strcmp(argv[0], "-c")) || (!std::strcmp(argv[0], "-contours"))) {
      opt.contour = true;
      argv++;
      argc--;
      std::cout << "-c ";
    }
    else if ((!std::strcmp(argv[0], "-s")) || (!std::strcmp(argv[0], "-shuffle"))) {
      opt.shuffle = true;
      argv++;
      argc--;
      std::cout << "-s ";
    }
    else if ((!std::strcmp(argv[0], "-b")) || (!std::strcmp(argv[0], "-binary"))) {
      opt.binary = true;
      argv++;
      argc--;
      std::cout << "-b ";
    }
    else if ((!std::strcmp(argv[0], "-x")) || (!std::strcmp(argv[0], "-xyz"))) {
      opt.xyz = true;
      argv++;
      argc--;
      std::cout << "-x ";
    }
    else if ((!std::strcmp(argv[0], "-nb")) || (!std::strcmp(argv[0], "-no_border"))) {
      opt.K = HUGE_VAL;
      opt.K_init = opt.K;
      argv++;
      argc--;
      std::cout << "-nb ";
    }
    else if ((!std::strcmp(argv[0], "-nh")) || (!std::strcmp(argv[0], "-no_header"))) {
      opt.no_header = true;
      argv++;
      argc--;
      std::cout << "-nh ";
    }
    else if ((!std::strcmp(argv[0], "-d")) || (!std::strcmp(argv[0], "-delta"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.delta) != 1) {
	std::cerr << "Argument for delta must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-d " << opt.delta << " "; 
    }  
    else if ((!std::strcmp(argv[0], "-a")) || (!std::strcmp(argv[0], "-area"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.area) != 1) {
	std::cerr << "Argument for area must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-a " << opt.area << " "; 
    }   
    else if ((!std::strcmp(argv[0], "-pe")) || (!std::strcmp(argv[0], "-perimeter"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.perimeter) != 1) {
	std::cerr << "Argument for perimeter must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-perimeter " << opt.perimeter << " "; 
    }  
    else if ((!std::strcmp(argv[0], "-ki"))){
      if ((CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K_init) != 1)||
	  (CGAL_CLIB_STD::sscanf(argv[2], "%lf", &opt.K) != 1)){
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 3;
      argc -= 3;
      std::cout << "-ki " << opt.K_init << " " << opt.K << " ";
    }  
    else if ((!std::strcmp(argv[0], "-rgb"))){
      if ((CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.red) != 1)||
	  (CGAL_CLIB_STD::sscanf(argv[2], "%lf", &opt.green) != 1) ||
	  (CGAL_CLIB_STD::sscanf(argv[3], "%lf", &opt.blue) != 1)){
	std::cerr << "Argument for rgb must be three numbers"
		  << std::endl;
      }
      argv += 4;
      argc -= 4;
      std::cout << "-rgb " << opt.red << " " << opt.green << " " << opt.blue << " " ;
    }
    else if ((!std::strcmp(argv[0], "-ks"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K_step) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-ks " << opt.K_step << " ";
    }  
    else if ((!std::strcmp(argv[0], "-k"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%lf", &opt.K) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      opt.K_init = opt.K;
      argv += 2;
      argc -= 2;
      std::cout << "-k " << opt.K_init << " ";
    } 
    else if ((!std::strcmp(argv[0], "-n")) || (!std::strcmp(argv[0], "-number_of_points"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.number_of_points) != 1) {
	std::cerr << "Argument for the number of points must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-n " << opt.number_of_points << " ";
    }  
    else if ((!std::strcmp(argv[0], "-m")) || (!std::strcmp(argv[0], "-max_of_connected_components"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.max_connected_comp) != 1) {
	std::cerr << "Argument for the number of connected components must be a number"
		  << std::endl;
      }
      /*
      if(opt.max_connected_comp < 1) {
	std::cerr << "Argument for the number of connected components must be a positive number"
		  << "It is set to 1" << std::endl;
	opt.max_connected_comp = 1;
      }
      */
      argv += 2;
      argc -= 2;
      std::cout << "-m " << opt.max_connected_comp << " ";
    }
    else if ((!std::strcmp(argv[0], "-p")) || (!std::strcmp(argv[0], "-post"))){
      if (CGAL_CLIB_STD::sscanf(argv[1], "%d", &opt.NB_BORDER_MAX) != 1) {
	std::cerr << "Argument for post process must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cout << "-p " << opt.NB_BORDER_MAX << " ";
    }
    else if ((!std::strcmp(argv[0], "-i")) || (!std::strcmp(argv[0], "-in"))) {
      std::strcpy(opt.finname, argv[1]);
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cout << "-i " << opt.finname << " ";
    }
    else if ((!std::strcmp(argv[0], "-s")) || (!std::strcmp(argv[0], "-sect_in"))) {
      std::strcpy(opt.finname, argv[1]);
      opt.Section_file = true;
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cout << "-s " << opt.finname << " ";
    }
    else if ((!std::strcmp(argv[0], "-o")) || (!std::strcmp(argv[0], "-out"))) {
      std::strcpy(opt.foutname, argv[1]);
      opt.file_output = true;
      argv += 2;
      argc -= 2;
      std::cout << "-o " << opt.foutname << " ";
    }
    else if ((!std::strcmp(argv[0], "-of")) || (!std::strcmp(argv[0], "-out_format"))) {
      if (!std::strcmp(argv[1], "wrl"))
	opt.out_format = 0;
      else if (!std::strcmp(argv[1], "oogl"))
	opt.out_format = 1;
      else if (!std::strcmp(argv[1], "medit"))
	opt.out_format = 2;
      else if (!std::strcmp(argv[1], "ply"))
	opt.out_format = 3;
      else if(!std::strcmp(argv[1], "iv"))
	opt.out_format = 4;
      else if(!std::strcmp(argv[1], "stl"))
	opt.out_format = 5;
      else if (!std::strcmp(argv[1], "all"))
	opt.out_format = -1;
      else if (!std::strcmp(argv[1], "none"))
	opt.out_format = -2;
      else
	std::cout << "unrecognized file format." << std::endl;
      opt.file_output = true;
      std::cout << "-of " << argv[1] << " ";
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

} // namespace NUAGE
