// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>

// typedefs for defining the adaptor
typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<K>               Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>                              DT;
typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<DT>                 AT;
typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<DT>         AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;

typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;

void print_endpoint(Halfedge_handle e, bool is_src) {
  std::cout << "\t";
  if ( is_src ) {
    if ( e->has_source() )  std::cout << e->source()->point() << std::endl;
    else  std::cout << "point at infinity" << std::endl;
  } else {
    if ( e->has_target() )  std::cout << e->target()->point() << std::endl;
    else  std::cout << "point at infinity" << std::endl;
  }
}

int main()
{
  std::ifstream ifs("data/data1.svd.cin");
  assert( ifs );

  VD vd;

  Site_2 t;
  while ( ifs >> t ) { vd.insert(t); }
  ifs.close();

  assert( vd.is_valid() );

  std::ifstream ifq("data/queries1.svd.cin");
  assert( ifq );

  Point_2 p;
  while ( ifq >> p ) {
    std::cout << "Query point (" << p.x() << "," << p.y()
              << ") lies on a Voronoi " << std::flush;

    Locate_result lr = vd.locate(p);
    if ( Vertex_handle* v = boost::get<Vertex_handle>(&lr) ) {
      std::cout << "vertex." << std::endl;
      std::cout << "The Voronoi vertex is:" << std::endl;
      std::cout << "\t" << (*v)->point() << std::endl;
    } else if ( Halfedge_handle* e = boost::get<Halfedge_handle>(&lr) ) {
      std::cout << "edge." << std::endl;
      std::cout << "The source and target vertices "
                << "of the Voronoi edge are:" << std::endl;
      print_endpoint(*e, true);
      print_endpoint(*e, false);
    } else if ( Face_handle* f = boost::get<Face_handle>(&lr) ) {
      std::cout << "face." << std::endl;
      std::cout << "The vertices of the Voronoi face are"
                << " (in counterclockwise order):" << std::endl;
      Ccb_halfedge_circulator ec_start = (*f)->ccb();
      Ccb_halfedge_circulator ec = ec_start;
      do {
        print_endpoint(ec, false);
      } while ( ++ec != ec_start );
    }
    std::cout << std::endl;
  }
  ifq.close();

  return 0;
}
