// point_location.C

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>

#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include <CGAL/Segment_Voronoi_diagram_Voronoi_traits_2.h>

typedef CGAL::Simple_cartesian<double>  K;

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_2<K>   Gt;
typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>        SVD;
typedef CGAL::Segment_Voronoi_diagram_Voronoi_traits_2<SVD>  VT;
typedef CGAL::Voronoi_diagram_adaptor_2<SVD,VT>              VD;

typedef VT::Assign                 Assign;
typedef VT::Point_locator::Object  Object;

typedef VD::Vertex_handle    Vertex_handle;
typedef VD::Face_handle      Face_handle;
typedef VD::Halfedge_handle  Halfedge_handle;



int main()
{
  std::ifstream ifs("data/data1.svd.cin");
  assert( ifs );

  SVD svd;

  std::cout << "Input sites:" << std::endl;
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

  Assign assign = vd.voronoi_traits().assign_object();

  while ( ifq >> p ) {
    std::cout << p << "\t --> \t" << std::flush;

    Object o = vd.locate(p);

    Vertex_handle      v;
    Face_handle        f;
    Halfedge_handle    e;

    if ( assign(v, o) ) {
      std::cout << "VERTEX";
    } else if ( assign(e, o) ) {
      std::cout << "EDGE";
    } else if ( assign(f, o) ) {
      std::cout << "FACE";
    }
    std::cout << std::endl;
  }

  return 0;
}
