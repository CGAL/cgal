#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <boost/shared_ptr.hpp>
#include <CGAL/gmp.h>
#include <CGAL/Timer.h>

template<class K>
void print_point ( CGAL::Point_2<K> const& p )
{
  //CGAL::exact(p);
  std::cout <<  p.x() << " " << p.y();
}

template<class K>
void print_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss )
{
  typedef CGAL::Straight_skeleton_2<K> Ss ;

  typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
  typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
  typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

  Halfedge_const_handle null_halfedge ;
  Vertex_const_handle   null_vertex ;

  std::cerr << "Straight skeleton with " << ss.size_of_vertices()
            << " vertices, " << ss.size_of_halfedges()
            << " halfedges and " << ss.size_of_faces()
            << " faces" << std::endl ;

  for ( Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i )
  {
    if ( !i->is_bisector() ) continue;
    std::cout << "2 " ;
    print_point(i->opposite()->vertex()->point()) ;
    std::cout << " 0 " ;
    print_point(i->vertex()->point());
    std::cout << " 0\n" ;
    //std::cout << " " << ( i->is_bisector() ? "bisector" : "contour" ) << std::endl;
  }
}



typedef CGAL::Exact_predicates_exact_constructions_kernel       EPEC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel       EPIC;
typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > SCLGQ;
typedef EPEC::Exact_kernel::FT                                EPEC_exact_ft;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt       EPEC_with_sqrt;



template< class Kernel, class FT>
void build_skeleton(const char* fname)
{
  typedef typename Kernel::Point_2                                         Point_2;
  typedef CGAL::Polygon_2<Kernel>                                 Polygon_2;
  typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>        SsBuilderTraits;
  typedef CGAL::Straight_skeleton_2<Kernel>                       Ss;
  typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss>   SsBuilder;

  Polygon_2 pgn;

  std::ifstream input(fname);

  FT x,y;
  while(input)
  {
    input >> x;
    if (!input) break;
    input >> y;
    if (!input) break;
    pgn.push_back( Point_2( typename Kernel::FT(x), typename Kernel::FT(y) ) );
  }
  input.close();

  std::cout << "Polygon has " << pgn.size() <<  " points\n";

  if(!pgn.is_counterclockwise_oriented()) {
      std::cerr << "Polygon is not CCW Oriented" << std::endl;
  }
  if(!pgn.is_simple()) {
      std::cerr << "Polygon is not simple" << std::endl;
  }

  CGAL::Timer time;
  time.start();
  SsBuilder ssb;
  ssb.enter_contour(pgn.vertices_begin(), pgn.vertices_end());
  boost::shared_ptr<Ss> straight_ske = ssb.construct_skeleton();
  time.stop();

  std::cout << "Time spent to build skeleton " << time.time() << "\n";

  if(!straight_ske->is_valid()) {
      std::cerr << "Straight skeleton is not valid" << std::endl;
  }

  std::cerr.precision(60);
  print_straight_skeleton(*straight_ske);

}




int main(int argc, char** argv) {
  if (argc!=2)
    std::cerr << "Usage: " << argv[0] << " polygon_points.txt\n";
  else
  {

    std::cout <<"Running with EPIC\n";
    build_skeleton<EPIC, double>(argv[1]);

    std::cout <<"Running with SCLGQ\n";
    build_skeleton<SCLGQ, SCLGQ::FT>(argv[1]);

    std::cout <<"Running with EPEC_with_sqrt\n";
    build_skeleton<EPEC_with_sqrt, EPEC_with_sqrt::FT>(argv[1]);

    std::cout <<"Running with EPEC\n";
    build_skeleton<EPEC, EPEC_exact_ft>(argv[1]);
  }

  return 0;
}
