#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K>             Gt;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<Gt>                Vbb;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Gt, Vbb>  Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2<Gt>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                    Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt, Tds>              Delaunay;
typedef Delaunay::Point                                                 Point;

//a functor that returns a std::pair<Point,unsigned>.
//the unsigned integer is incremented at each call to
//operator()
struct Auto_count
  : public CGAL::cpp98::unary_function<const Point&, std::pair<Point, unsigned> >
{
  mutable unsigned i;
  Auto_count() : i(0) {}
  std::pair<Point, unsigned> operator()(const Point& p) const
  {
    return std::make_pair(p, i++);
  }
};

int main()
{
  std::vector<Point> points;
  points.push_back( Point(0.0, 0.0) );
  points.push_back( Point(0.1, 0.0) );
  points.push_back( Point(0.0, 0.1) );
  points.push_back( Point(0.1, 0.4) );
  points.push_back( Point(0.2, 0.2) );
  points.push_back( Point(0.4, 0.0) );

  Delaunay T;
  T.insert( boost::make_transform_iterator(points.begin(), Auto_count()),
            boost::make_transform_iterator(points.end(),  Auto_count() )  );

  CGAL_assertion( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    if( points[ vit->info() ] != vit->point() )
      {
        std::cerr << "Error different info" << std::endl;
        exit(EXIT_FAILURE);
      }
  }
  std::cout << "OK" << std::endl;

  return 0;
}
