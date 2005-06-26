// examples/Voronoi_diagram_adaptor_2/point_location.C

#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>
#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include <CGAL/Segment_Voronoi_diagram_Voronoi_traits_2.h>

// typedefs for defining the adaptor
typedef CGAL::Simple_cartesian<double>                       K;
typedef CGAL::Segment_Voronoi_diagram_filtered_traits_2<K>   Gt;
typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>        SVD;
typedef CGAL::Segment_Voronoi_diagram_Voronoi_traits_2<SVD>  VT;
typedef CGAL::Voronoi_diagram_adaptor_2<SVD,VT>              VD;

// typedef for the result type of the point location
typedef VD::Locate_result      Locate_result;


int main()
{
  std::ifstream ifs("data/data1.svd.cin");
  assert( ifs );

  std::cout << "Input sites:" << std::endl;

  SVD svd;
  SVD::Site_2 t;
  while ( ifs >> t ) {
    std::cout << t << std::endl;
    svd.insert(t);
  }
  std::cout << std::endl;
  ifs.close();

  assert( svd.is_valid() );

  if ( svd.dimension() < 1 ) {
    std::cout << "The number of generators is less than 2." << std::endl;
    std::cout << "Cannot do point location." << std::endl;
    return 0;
  }

  VD vd(svd);
  
  std::ifstream ifq("data/queries1.svd.cin");
  assert( ifq );

  std::cout << "Query sites and location feature:" << std::endl;

  SVD::Point_2 p;
  while ( ifq >> p ) {
    std::cout << p << "\t --> \t" << std::flush;

    Locate_result lr = vd.locate(p);
    if ( lr.is_vertex() ) {
      std::cout << "VERTEX";
    } else if ( lr.is_edge() ) {
      std::cout << "EDGE";
    } else if ( lr.is_face() ) {
      std::cout << "FACE";
    }
    std::cout << std::endl;
  }

  return 0;
}
