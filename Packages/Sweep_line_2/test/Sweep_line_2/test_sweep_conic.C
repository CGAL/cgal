#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Arr_conic_traits_2.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Bench_parse_args.h>
#define USE_CONIC_TRAITS
#include <../../../Arrangement/bench/Arrangement/Conic_reader.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>


typedef leda_real                         NT;
typedef CGAL::Cartesian<NT>               Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>  Traits;

typedef Traits::X_monotone_curve_2        X_monotone_curve_2;
typedef Traits::Point_2                   Point_2;

typedef CGAL::Pm_default_dcel<Traits>                    Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                  Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>        Arr_2;

typedef std::list<X_monotone_curve_2>     CurveContainer;
typedef CurveContainer::iterator CurveContainerIter;

typedef Conic_reader<Traits> MyReader;
typedef std::list<Point_2> PointList;
typedef PointList::iterator PointListIter;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
CGAL::Pm_naive_point_location<Arr_2::Planar_map> pl; 

static Arr_2               arr(&pl);



//---------------------------------------------------------------------------
// The main:
//
int main (int argc, char** argv)
{
  bool verbose = false;

  // Define a test objects to read the conic arcs from it.
  if (argc<2) 
  {
    std::cerr << "Usage: Conic_traits_test <filename>" << std::endl;
    exit(1);
  }

  CGAL::Bbox_2 bbox;
  CurveContainer curves;

  MyReader reader;
  reader.ReadData(argv[1], curves, CGAL::Bench_parse_args::FORMAT_INT, bbox);

  // run the sweep
  std::list<X_monotone_curve_2> mylist;

  CGAL::Sweep_line_2<CurveContainerIter, Traits> sl;
  sl.get_subcurves(curves.begin(), curves.end(), 
		   std::back_inserter(mylist), false);
  
  
  PointList point_list_with_ends;
  CGAL::Sweep_line_2<CurveContainerIter, Traits> sl1;
  sl1.get_intersection_points(curves.begin(), curves.end(), 
			      std::back_inserter(point_list_with_ends));
  int point_count_with_ends_calculated = point_list_with_ends.size();
  
  // generate the string for the output
  std::stringstream out1;
  for ( std::list<X_monotone_curve_2>::iterator iter = mylist.begin() ;
	iter != mylist.end() ; ++iter )
  {
    out1 << *iter << "\n";
  }
  
  // read the output from the file
  std::stringstream out2;
  char buf[1024];
  int count = 0;
  
  std::ifstream in_file(argv[1]);
  in_file >> count;
  in_file.getline(buf, 1024); // to get rid of the new line
  for ( int i = 0 ; i < count ; i++ ) {
    in_file.getline(buf, 1024);
  }
  in_file >> count;
  in_file.getline(buf, 1024); // to get rid of the new line
  for (int i = 0; i < count; i++) {
    in_file.getline(buf, 1024);
    out2 << buf << "\n";
  }
  int point_count_with_ends_from_file = 0;
  in_file >> point_count_with_ends_from_file;
  in_file.close();
  
  if ( verbose )
  {
    std::cout << "Result: \n" << mylist.size() << "\n";
    for ( std::list<X_monotone_curve_2>::iterator i = mylist.begin() ;
	  i != mylist.end() ; ++i )
    {
      std::cout << *i << "\n";
    }
  }
  
  std::string calculated = out1.str();
  std::string infile = out2.str();
  
  if ( infile == calculated ) {
    if ( point_count_with_ends_from_file != 
	 point_count_with_ends_calculated ) {
      std::cout << "number of intersection points (with ends):" 
		<< point_count_with_ends_calculated << ". Should be " 
		<< point_count_with_ends_from_file << "\n";
      std::cout << argv[1] << " Error\n";
      return -1;
    }  else {
      std::cout << argv[1] << " OK!\n";
    }
  } else {
    std::cout << argv[1] << " Error\n";
    std::cout << "\ncalculated:\n";
    std::cout << calculated << std::endl;
    std::cout << "\nin file:\n";
    std::cout << infile << std::endl;
    std::cout << "--"  << std::endl;
    return -1;
  }
  
  return 0;  
}


