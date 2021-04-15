#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Vertex_handle                                     Vertex_handle;

int main()
{
  std::vector< std::pair<Point,unsigned> > points;
  points.push_back( std::make_pair(Point(0,0),0)   );
  points.push_back( std::make_pair(Point(1,0),1)   );
  points.push_back( std::make_pair(Point(0,1),2)   );
  points.push_back( std::make_pair(Point(14,4),3)  );
  points.push_back( std::make_pair(Point(2,2),4)   );
  points.push_back( std::make_pair(Point(-4,0),5)  );


  Delaunay T;
  T.insert( points.begin(),points.end() );

  CGAL_assertion( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (Vertex_handle v : T.finite_vertex_handles())
    if( points[ v->info() ].first != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
