#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point                                             Point;

int main()
{
  std::vector< std::pair<Point,unsigned> > points;
  points.push_back( std::make_pair(Point(0,0,0),0) );
  points.push_back( std::make_pair(Point(1,0,0),1) );
  points.push_back( std::make_pair(Point(0,1,0),2) );
  points.push_back( std::make_pair(Point(0,0,1),3) );
  points.push_back( std::make_pair(Point(2,2,2),4) );
  points.push_back( std::make_pair(Point(-1,0,1),5) );


  Delaunay T( points.begin(),points.end() );

  CGAL_assertion( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  for (Delaunay::Vertex_handle v : T.finite_vertex_handles())
    if( points[ v->info() ].first != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
