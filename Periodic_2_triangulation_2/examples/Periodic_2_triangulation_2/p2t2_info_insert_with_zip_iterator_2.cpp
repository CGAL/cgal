#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Periodic_2_triangulation_traits_2<K>                  Gt;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Gt>   Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2<Gt>              Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                 Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt, Tds>          Delaunay;
typedef Delaunay::Point                                             Point;

int main()
{

  std::vector<unsigned> indices;
  indices.push_back(0);
  indices.push_back(1);
  indices.push_back(2);
  indices.push_back(3);
  indices.push_back(4);
  indices.push_back(5);

  std::vector<Point> points;
  points.push_back( Point(0.0, 0.0) );
  points.push_back( Point(0.1, 0.0) );
  points.push_back( Point(0.0, 0.1) );
  points.push_back( Point(0.1, 0.4) );
  points.push_back( Point(0.2, 0.2) );
  points.push_back( Point(0.4, 0.0) );

  Delaunay T;
  T.insert( boost::make_zip_iterator(boost::make_tuple( points.begin(), indices.begin() )),
            boost::make_zip_iterator(boost::make_tuple( points.end(), indices.end() ) )  );

  CGAL_assertion( T.number_of_vertices() == 6 );


  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    if( points[ vit->info() ] != vit->point() )
      {
        std::cerr << "Error different info" << std::endl;
        exit(EXIT_FAILURE);
      }
  std::cout << "OK" << std::endl;

  return 0;
}
