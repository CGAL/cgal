#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;

typedef CGAL::Regular_triangulation_vertex_base_2<K>                   Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K,Vbase> Vb;
typedef CGAL::Regular_triangulation_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
typedef CGAL::Regular_triangulation_2<K, Tds>                          Regular;
typedef K::Point_2                                                     Point;
typedef K::Weighted_point_2                                            Wpoint;
typedef Regular::Vertex_handle                                         Vertex_handle;

int main()
{
  std::vector< std::pair<Wpoint,unsigned> > points;
  points.push_back( std::make_pair(Wpoint(Point(0,0),2),0) );
  points.push_back( std::make_pair(Wpoint(Point(1,0),2),1) );
  points.push_back( std::make_pair(Wpoint(Point(0,1),2),2) );
  points.push_back( std::make_pair(Wpoint(Point(-4,54),2),3) );
  points.push_back( std::make_pair(Wpoint(Point(2,2),2),4) );
  points.push_back( std::make_pair(Wpoint(Point(-1,0),2),5) );

  Regular rt;
  rt.insert( points.begin(),points.end() );

  CGAL_assertion( rt.number_of_vertices() == 6 );

  // check that the info was correctly set.
  for (Vertex_handle v : rt.finite_vertex_handles())
    if( points[ v->info() ].first != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
