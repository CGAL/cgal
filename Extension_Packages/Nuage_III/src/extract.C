//#define DEBUG
//#define VERBOSE
//#define FACET_NUMBER
//=====================================================================
// selection de facettes dans Delaunay....
//=====================================================================
#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>
#include <cassert>
#include <vector>
#include <list>


// Kernel
#include <CGAL/Simple_cartesian.h>

#if 0
#include <CGAL/Filtered_kernel.h>
#else
#include <CGAL/Static_filters.h>
#endif

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/Timer.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <NUAGE/Local_selection_vertex_base_3.h>

#include <CGAL/Triangulation_cell_base_3.h>
#include <NUAGE/Local_selection_cell_base_3.h>
#include <NUAGE/Extract_surface.h>

#include "Parse.C"

//=====================================================================
//=====================================================================

typedef double coord_type;
typedef double NT;

#if 0
struct Kernel : public CGAL::Filtered_kernel<CGAL::Simple_cartesian<NT> > {};
#else
struct Kernel : public CGAL::Static_filters<CGAL::Simple_cartesian<NT> > {};
#endif


typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_2  Point_2;

typedef CGAL::Triangulation_vertex_base_3<Kernel> TVb;
typedef Local_selection_vertex_base_3<TVb> LVb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<LVb> HVb;

typedef CGAL::Triangulation_cell_base_3<Kernel> Cb;
typedef Local_selection_cell_base_3<Cb> LCb;


typedef CGAL::Triangulation_data_structure_3<HVb,LCb> Tds;


typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Delaunay_Triangulation_3;
typedef CGAL::Triangulation_hierarchy_3<Delaunay_Triangulation_3> Triangulation_3;

typedef Extract_surface<Triangulation_3,Kernel> Surface;


CGAL::Timer t1;

#include <NUAGE/iofile_manipulator.h>

void
file_output(char* foutput, std::vector<Point>& points)
{
  std::ofstream os(foutput, std::ios::out);

  if(os.fail())
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cerr << "+++file for output : " << foutput << std::endl;

  os.clear();

  CGAL::set_ascii_mode(os);

  for(std::vector<Point>::iterator it = points.begin();
      it != points.end(); it++)
    os << *it;
}

//---------------------------------------------------------------------

bool
file_input(const Options& opt, std::vector<Point>& points)
{
  const char* finput = opt.finname;
  int number_of_points =  opt.number_of_points;
  bool xyz = opt.xyz;
  
  std::ios::openmode mode = (opt.binary) ? std::ios::binary : std::ios::in;
  std::ifstream is(finput, mode);

  if(opt.binary){
    CGAL::set_binary_mode(is);
  }
  if(is.fail())
    {
      std::cerr << "+++unable to open file for input" << std::endl;
      CGAL_CLIB_STD::exit(0);
      return false;
    }
  else
    std::cout << "Input from file : " << finput << std::endl;

// pour selectionner le mode de lecture souhaite...
//   is.setf(std::ifstream::scientific);
//   is.setf(std::ifstream::showpos);
//   is.setf(std::ifstream::uppercase);



  int n;
  if(! xyz){
    is >> n;
    std::cout << "   reading " << n << " points" << std::endl;
    points.reserve(n);
    CGAL::copy_n(std::istream_iterator<Point>(is), n, std::back_inserter(points));
  } else {
    // we do not know beforehand how many points we will read
    std::istream_iterator<Point> it(is), eof;
    while(it!= eof){
      points.push_back(*it);
      it++;
    }
    n = points.size();
  }

  if(opt.shuffle) {
    std::cout << "Random shuffling" << std::endl;
    std::random_shuffle(points.begin(), points.end());
  }


  if ( (number_of_points > 0 ) && (number_of_points < n ))
    {
      points.erase(points.begin()+number_of_points, points.begin()+n);

      std::cout << std::endl 
		<< "   and randomize a sub-sample of " << number_of_points 
		<< " points." <<
	std::endl << std::endl;
    }
  return true;
}

//---------------------------------------------------------------------

bool
section_file_input(char* finput, const int& number_of_points, std::vector<Point>& points)
{
  std::ifstream is(finput, std::ios::in);

  if(is.fail())
    {
      std::cerr << "+++unable to open file for input" << std::endl;
      CGAL_CLIB_STD::exit(0);
      return false;
    }
  else
    std::cout << ">> input from section-file : " << finput << std::endl;

// pour selectionner le mode de lecture souhaite...
//   is.setf(std::ifstream::scientific);
//   is.setf(std::ifstream::showpos);
//   is.setf(std::ifstream::uppercase);

  CGAL::set_ascii_mode(is);
      
  int i(0), N, points_num(0);
  char c;
  is >> c; // c == "S"
  is >> N;
  coord_type h;
  int n;

  do{
    is >> c; // c == "v"
    is >> n;
    is >> c; // c == "z"
    is >> h;
    is >> c; // c == "{"
    
    Point_2 p;

    for(; n > 0; n--)
      {
	is >> c;
	if (c == '}') 
	  {
	    is >> c; // c == "{"
	    n++;
	  }
	else
	  {
	    is.putback(c);
	    is >> p;
	    points_num++;
	    points.push_back(Point (p.x(), p.y(), h));
	  }
      }

    is >> c; // c == "}"
    i++;
  } while (i < N); 

  std::cout << "Reading " << points_num << " points";

  std::random_shuffle(points.begin(), points.end());

  if ( (number_of_points > 0 ) && (number_of_points < points_num ))
    {
      points.erase(points.begin()+number_of_points, points.begin()+points_num);

      std::cout << std::endl 
		<< "Randomize a sub-sample of " << number_of_points 
		<< " points." <<
	std::endl << std::endl;
    }
  return true;
}

void
construct_delaunay(const std::vector<Point> &V_p,
		   Triangulation_3& T)
{
  std::cout << "Compute Delaunay Tetrahedrization " << std::endl; 
  t1.start();
  {
    for(std::vector<Point>::const_iterator v_it = V_p.begin();
	v_it != V_p.end(); ++v_it)
      {
		  Triangulation_3::Vertex_handle vh = T.insert(*v_it);
      }
  }
  t1.stop();
  std::cout << "   Inserted " << T.number_of_vertices() << " points, "
	    <<  T.number_of_cells() << " cells computed in "
	    << t1.time() << " sec." << std::endl;
  if (T.dimension() < 3)
    {
      std::cout << "-- 2D sample of points ???"  << std::endl;
      CGAL_CLIB_STD::exit(0);
    }
  t1.reset();
}




//___________________________________________
int main(int argc,  char* argv[])
{
  CGAL::Timer t2;

  t2.start();
  //parse command line
  Options opt;
  std::cout << "Option line for this execution is :" << std::endl;
  if (!parse(argc, argv, opt))
    CGAL_CLIB_STD::exit(0);
  std::cout << std::endl << std::endl;
  
  Triangulation_3 T;

  std::vector<Point> points;

  if (!opt.Section_file)
    file_input(opt, points);
  else
    section_file_input(opt.finname, opt.number_of_points, points);

  std::cout << "Time for reading "  << t2.time() << " sec." << std::endl;
  construct_delaunay(points, T);
  
  points.clear();

  bool re_init = false;
  int number_of_connected_comp = 0;
  double  postprocessing_time = 0;

  Surface S(T,opt.DELTA);
  do 
    {
      number_of_connected_comp++;

      if (S.init(re_init))
	{
	  re_init = false;
	  std::cout << "Growing connected component " << number_of_connected_comp << std::endl;
	  S.extend(opt.K_init, opt.K_step, opt.K);

	  std::cout << "Number of facets: " << S.number_of_facets() << std::endl;
	  
	  if ((S.number_of_facets() > T.number_of_vertices())&&
	      (opt.NB_BORDER_MAX > 0))
	    // en principe 2*nb_sommets = nb_facettes: y a encore de la marge!!!
	    {
	      t1.reset();
	      while(S.postprocessing(opt.NB_BORDER_MAX)){
		re_init = true;
		S.extend(opt.K_init, opt.K_step, opt.K);
	      }
	    }

	    postprocessing_time += t1.time();

	} 
    }while(re_init &&
	   ((number_of_connected_comp < opt.max_connected_comp)||
	    (opt.max_connected_comp < 0)));

 
  std::cout << "Total time: " << t2.time() << " sec." << std::endl; 
  write_in_file_selected_facets(opt.foutname, S, opt.contour, opt.out_format);


  std::cout << "   "  << S.number_of_outliers()
	    << " outliers." << std::endl; 
  std::cout << "   Reconstructed surface: " << S.number_of_facets() << 
    " facets, " << S.number_of_vertices() << " vertices." << std::endl;
  std::cout << "   "  << S.number_of_border_edges() << 
    " border edges." << std::endl;
  std::cout << "   number of connected components <= " 
	    << std::max(1, number_of_connected_comp-1)
	    << std::endl << std::endl;

  return 0;
}


