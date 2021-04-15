#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;

typedef CGAL::Regular_triangulation_vertex_base_3<K>                    Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K, Vbase> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K>                      Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>                     Tds;
typedef CGAL::Regular_triangulation_3<K, Tds>                           Regular;
typedef K::Point_3                                                      Point;
typedef K::Weighted_point_3                                        Wpoint;

int main()
{
  std::vector< std::pair<Wpoint,unsigned> > points;
  points.push_back( std::make_pair(Wpoint(Point(0,0,0),2),0) );
  points.push_back( std::make_pair(Wpoint(Point(1,0,0),2),1) );
  points.push_back( std::make_pair(Wpoint(Point(0,1,0),2),2) );
  points.push_back( std::make_pair(Wpoint(Point(0,0,1),2),3) );
  points.push_back( std::make_pair(Wpoint(Point(2,2,2),2),4) );
  points.push_back( std::make_pair(Wpoint(Point(-1,0,1),2),5) );

  Regular rt( points.begin(),points.end() );

  CGAL_assertion( rt.number_of_vertices() == 6 );

  // check that the info was correctly set.
  for (Regular::Vertex_handle v : rt.finite_vertex_handles())
    if( points[ v->info() ].first != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
