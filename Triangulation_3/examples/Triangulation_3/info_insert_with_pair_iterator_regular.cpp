#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                     K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>                       Traits;

typedef CGAL::Exact_predicates_inexact_constructions_kernel                     K;
typedef CGAL::Triangulation_vertex_base_3<Traits>                               Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, Traits,Vbase>     Vb;
typedef CGAL::Regular_triangulation_cell_base_3<Traits>                         Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>                             Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds>                              Regular;
typedef K::Point_3                                                              Point;
typedef Traits::Weighted_point_3                                                Wpoint;

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
  Regular::Finite_vertices_iterator vit;
  for (vit = rt.finite_vertices_begin(); vit != rt.finite_vertices_end(); ++vit)
    if( points[ vit->info() ].first != vit->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
