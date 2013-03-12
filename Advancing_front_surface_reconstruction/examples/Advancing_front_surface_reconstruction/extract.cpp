#define NOLAZY
#define BLIND


#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>
#include <cassert>
#include <vector>
#include <list>

#include <boost/iterator/transform_iterator.hpp>

// Kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/AFSR_vertex_base_with_id_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/AFSR_cell_base_3.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/AFSR_options.h>
#include <CGAL/IO/Advancing_front_surface_reconstruction.h>





typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_3  Point;

typedef CGAL::AFSR_vertex_base_with_id_3<Kernel> LVb;

typedef CGAL::Triangulation_cell_base_3<Kernel> Cb;
typedef CGAL::AFSR_cell_base_3<Cb> LCb;

typedef CGAL::Triangulation_data_structure_3<LVb,LCb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

typedef Triangulation_3::Vertex_handle Vertex_handle;

typedef CGAL::Advancing_front_surface_reconstruction<Kernel,Triangulation_3> Surface;
typedef CGAL::AFSR_options Options;


//---------------------------------------------------------------------

struct Auto_count : public std::unary_function<const Point&,std::pair<Point,int> >{
  mutable int i;
  Auto_count() : i(0){}
  std::pair<Point,int> operator()(const Point& p) const {
    return std::make_pair(p,i++);
  }
};


bool
file_input(const Options& opt, std::vector<Point>& points)
{
  const char* finput = opt.finname;
  bool xyz = opt.xyz;
  
  std::ios::openmode mode = (opt.binary) ? std::ios::binary : std::ios::in;
  std::ifstream is(finput, mode);

  if(opt.binary){
    CGAL::set_binary_mode(is);
  }
  if(is.fail())
    {
      std::cerr << "+++unable to open file for input" << std::endl;
      exit(0);
      return false;
    }
  else
    std::cerr << "Input from file : " << finput << std::endl;

  std::size_t n;
  if(! xyz){
    is >> n;
    std::cerr << "   reading " << n << " points" << std::endl;
    points.reserve(n);
     CGAL::cpp11::copy_n(std::istream_iterator<Point>(is), n, std::back_inserter(points));
  } else {
    // we do not know beforehand how many points we will read
    std::istream_iterator<Point> it(is), eof;
    char ignore[256];
    while(it!= eof){
      points.push_back(*it);
      is.getline(ignore,256); // ignore what comes after 3 doubles in a line
      it++;
    }
    n = points.size();
  }

  return true;
}



void usage(char* program)
{
  std::cerr << std::endl << "NAME     " << std::endl
	    << program << "     - surface extension -" << std::endl << std::endl;

  std::cerr << std::endl << "OPTIONS" << std::endl
	    << "     -xyz             : input data in xyz format" << std::endl
	    << "     -no_border -nb   : " << std::endl
	    << "     -in fname        : reads points from file ./fname" << std::endl
	    << "     -out fname       : writes points to file ./fname" << std::endl
	    << "     -out_format -of  : choose file format for output (iv, wrl, off, medit," << std::endl
	    << "                                                     ply, stl, all, none)" << std::endl
            << "     -rgb r g b       : color of the surface" << std::endl
            << "     -no_header       : The Vrml header and footer are not written" << std::endl	    
	    << "     -area  a         : No faces larger than area * average_area" << std::endl
	    << "     -perimeter  p    : No faces larger than perimeter * average_perimeter" << std::endl   
	    << "     -abs_area  a     : No faces larger than abs_area" << std::endl
	    << "     -abs_perimeter p : No faces with perimeter longer than abs_perimeter" << std::endl
            << "\n     Options for internal use" << std::endl

	    << "     -sect_in fname   : reads points from sections file ./fname" << std::endl
	    << "     -binary          : binary I/O" << std::endl
	    << "     -delta x         : set the delta constant" << std::endl
	    << "     -ki x y          : set the K interval (default : [1.1  5])" << std::endl
	    << "     -ks x            : set the K step (default : .1)" << std::endl
	    << "     -k  x            : set the K constant (only one pass)" << std::endl
	    << "     -Delaunay        : display the underlying Delaunay triangulation" << std::endl
	    << "     -max_of_connected_components x     : set the max of connected components" << std::endl
	    << "                                          (default : non-active)" << std::endl
            << "     -post x          : set a number for the post process" << std::endl
	    << "     -contours        : display contours" << std::endl;
}



bool
parse(int argc, char* argv[], Options &opt)
{
  std::strcpy(opt.program, argv[0]);
  --argc;
  argv++;
  if(argc == 0)
    std::cerr << "nothing ???" << std::endl;
  
  while ((argc > 0) && (argv[0][0] == '-')){
    if ((!std::strcmp(argv[0], "-D")) || (!std::strcmp(argv[0], "-Delaunay"))) {
      opt.Delaunay = true;
      argv++;
      argc--;
      std::cerr << "-D ";
    }
    else if ((!std::strcmp(argv[0], "-c")) || (!std::strcmp(argv[0], "-contours"))) {
      opt.contour = true;
      argv++;
      argc--;
      std::cerr << "-c ";
    }
    else if ((!std::strcmp(argv[0], "-b")) || (!std::strcmp(argv[0], "-binary"))) {
      opt.binary = true;
      argv++;
      argc--;
      std::cerr << "-b ";
    }
    else if ((!std::strcmp(argv[0], "-x")) || (!std::strcmp(argv[0], "-xyz"))) {
      opt.xyz = true;
      argv++;
      argc--;
      std::cerr << "-x ";
    }
    else if ((!std::strcmp(argv[0], "-nb")) || (!std::strcmp(argv[0], "-no_border"))) {
      opt.K = HUGE_VAL;
      opt.K_init = opt.K;
      argv++;
      argc--;
      std::cerr << "-nb ";
    }
    else if ((!std::strcmp(argv[0], "-nh")) || (!std::strcmp(argv[0], "-no_header"))) {
      opt.no_header = true;
      argv++;
      argc--;
      std::cerr << "-nh ";
    }
    else if ((!std::strcmp(argv[0], "-d")) || (!std::strcmp(argv[0], "-delta"))){
      if (sscanf(argv[1], "%lf", &opt.delta) != 1) {
	std::cerr << "Argument for delta must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-d " << opt.delta << " "; 
    }  
    else if ((!std::strcmp(argv[0], "-a")) || (!std::strcmp(argv[0], "-area"))){
      if (sscanf(argv[1], "%lf", &opt.area) != 1) {
	std::cerr << "Argument for area must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-a " << opt.area << " "; 
    }   
    else if ((!std::strcmp(argv[0], "-pe")) || (!std::strcmp(argv[0], "-perimeter"))){
      if (sscanf(argv[1], "%lf", &opt.perimeter) != 1) {
	std::cerr << "Argument for perimeter must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-perimeter " << opt.perimeter << " "; 
    }    
    else if ((!std::strcmp(argv[0], "-aa")) || (!std::strcmp(argv[0], "-abs_area"))){
      if (sscanf(argv[1], "%lf", &opt.abs_area) != 1) {
	std::cerr << "Argument for abs_area must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-abs_area " << opt.abs_area << " "; 
    }     
    else if ((!std::strcmp(argv[0], "-ae")) || (!std::strcmp(argv[0], "-abs_perimeter"))){
      if (sscanf(argv[1], "%lf", &opt.abs_perimeter) != 1) {
	std::cerr << "Argument for abs_perimeter must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-abs_perimeter " << opt.abs_perimeter << " "; 
    }  
    else if ((!std::strcmp(argv[0], "-ki"))){
      if ((sscanf(argv[1], "%lf", &opt.K_init) != 1)||
	  (sscanf(argv[2], "%lf", &opt.K) != 1)){
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 3;
      argc -= 3;
      std::cerr << "-ki " << opt.K_init << " " << opt.K << " ";
    }  
    else if ((!std::strcmp(argv[0], "-rgb"))){
      if ((sscanf(argv[1], "%lf", &opt.red) != 1)||
	  (sscanf(argv[2], "%lf", &opt.green) != 1) ||
	  (sscanf(argv[3], "%lf", &opt.blue) != 1)){
	std::cerr << "Argument for rgb must be three numbers"
		  << std::endl;
      }
      argv += 4;
      argc -= 4;
      std::cerr << "-rgb " << opt.red << " " << opt.green << " " << opt.blue << " " ;
    }
    else if ((!std::strcmp(argv[0], "-ks"))){
      if (sscanf(argv[1], "%lf", &opt.K_step) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-ks " << opt.K_step << " ";
    }  
    else if ((!std::strcmp(argv[0], "-k"))){
      if (sscanf(argv[1], "%lf", &opt.K) != 1) {
	std::cerr << "Argument for K must be a number"
		  << std::endl;
      }
      opt.K_init = opt.K;
      argv += 2;
      argc -= 2;
      std::cerr << "-k " << opt.K_init << " ";
    } 
    else if ((!std::strcmp(argv[0], "-m")) || (!std::strcmp(argv[0], "-max_of_connected_components"))){
      if (sscanf(argv[1], "%d", &opt.max_connected_comp) != 1) {
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
      std::cerr << "-m " << opt.max_connected_comp << " ";
    }
    else if ((!std::strcmp(argv[0], "-p")) || (!std::strcmp(argv[0], "-post"))){
      if (sscanf(argv[1], "%d", &opt.NB_BORDER_MAX) != 1) {
	std::cerr << "Argument for post process must be a number"
		  << std::endl;
      }
      argv += 2;
      argc -= 2;
      std::cerr << "-p " << opt.NB_BORDER_MAX << " ";
    }
    else if ((!std::strcmp(argv[0], "-i")) || (!std::strcmp(argv[0], "-in"))) {
      std::strcpy(opt.finname, argv[1]);
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cerr << "-i " << opt.finname << " ";
    }
    else if ((!std::strcmp(argv[0], "-s")) || (!std::strcmp(argv[0], "-sect_in"))) {
      std::strcpy(opt.finname, argv[1]);
      opt.Section_file = true;
      opt.file_input = true;
      argv += 2;
      argc -= 2;
      std::cerr << "-s " << opt.finname << " ";
    }
    else if ((!std::strcmp(argv[0], "-o")) || (!std::strcmp(argv[0], "-out"))) {
      std::strcpy(opt.foutname, argv[1]);
      opt.file_output = true;
      argv += 2;
      argc -= 2;
      std::cerr << "-o " << opt.foutname << " ";
    }
    else if ((!std::strcmp(argv[0], "-of")) || (!std::strcmp(argv[0], "-out_format"))) {
      if (!std::strcmp(argv[1], "wrl"))
	opt.out_format = 0;
      else if (!std::strcmp(argv[1], "off"))
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
	std::cerr << "unrecognized file format." << std::endl;
      opt.file_output = true;
      std::cerr << "-of " << argv[1] << " ";
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

template <class PointIterator, class TripleOutputIterator>
void reconstruction_test(PointIterator point_begin, PointIterator
                         point_end, TripleOutputIterator out,  bool filter_input_points=false,
                         double perimeter=0) 
{
  Options opt;
  opt.abs_perimeter = perimeter;
  std::cerr << "Compute Delaunay Tetrahedrization " << std::endl; 
  CGAL::Timer t1;
  t1.start();
  
  Triangulation_3 dt( boost::make_transform_iterator(point_begin, Auto_count()),
                      boost::make_transform_iterator(point_end, Auto_count() )  );
  t1.stop();
  std::cerr << "   Inserted " << dt.number_of_vertices() << " points, "
	    <<  dt.number_of_cells() << " cells computed in "
	    << t1.time() << " sec." << std::endl;

  Surface S(dt, opt);

  write_triple_indices(out, S);
}


//___________________________________________
int main(int argc,  char* argv[])
{
  CGAL::Timer timer, total;
  total.start();
  timer.start();
  //parse command line
  Options opt;
  std::cerr << "Option line for this execution is :" << std::endl;
  if (!parse(argc, argv, opt))
    exit(0);
  std::cerr << std::endl << std::endl;
  
  std::vector<Point> points;

  file_input(opt, points);
 
  std::cerr << "Time for reading "  << timer.time() << " sec." << std::endl;
  std::vector<CGAL::Triple<int,int,int> > triples;
  reconstruction_test(points.begin(), points.end(), std::back_inserter(triples));
  std::cerr << triples.size() << std::endl;
#if 0
  std::cerr << "Compute Delaunay Tetrahedrization " << std::endl; 
  CGAL::Timer t1;
  t1.start();
  
  Triangulation_3 dt( boost::make_transform_iterator(points.begin(),Auto_count()),
                      boost::make_transform_iterator(points.end(),  Auto_count() )  );
  t1.stop();
  std::cerr << "   Inserted " << dt.number_of_vertices() << " points, "
	    <<  dt.number_of_cells() << " cells computed in "
	    << t1.time() << " sec." << std::endl;
 
  if (dt.dimension() < 3) {
    std::cerr << "-- 2D sample of points ???"  << std::endl;
    exit(0);
  }
  
  points.clear();

  
  Surface S(dt, opt);


  std::cerr << "Total time: " << timer.time() << " sec." << std::endl; 
  //  write_to_file_vrml2(opt.foutname, S, opt.contour, opt.red, opt.green, opt.blue, opt.no_header);
  //  write_to_file(opt.foutname, S, opt.contour, opt.out_format, opt.red, opt.green, opt.blue, opt.no_header);
  std::vector<CGAL::Triple<int,int,int> > triples;
  write_triple_indices(std::back_inserter(triples), S);
  std::cerr << triples.size() << std::endl;

  std::cerr << "   "  << S.number_of_outliers()
	    << " outliers." << std::endl; 
  std::cerr << "   Reconstructed surface: " << S.number_of_facets() << 
    " facets, " << S.number_of_vertices() << " vertices." << std::endl;
  std::cerr << "   "  << S.number_of_border_edges() << 
    " border edges." << std::endl;
  std::cerr << "   number of connected components <= " 
	    << (std::max)(1, S.number_of_connected_components()-1)
	    << std::endl << std::endl;
  
#endif
  total.stop();
  std::cerr << "Total = " << total.time() << " sec." << std::endl;
  return 0;
}




