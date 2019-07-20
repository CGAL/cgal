#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <boost/iterator/zip_iterator.hpp>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
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
  points.push_back(Point(0,0,0));
  points.push_back(Point(1,0,0));
  points.push_back(Point(0,1,0));
  points.push_back(Point(0,0,1));
  points.push_back(Point(2,2,2));
  points.push_back(Point(-1,0,1));

  Delaunay T( boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )  );

  CGAL_assertion( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  for (Delaunay::Vertex_handle v : T.finite_vertex_handles() )
    if( points[ v->info() ] != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }

  return 0;
}
