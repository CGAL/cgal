#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Triangulation_cell_base_3<K>                          Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
typedef CGAL::Triangulation_3<K, Tds>                               Triangulation;
typedef Triangulation::Point                                        Point;

//a functor that returns a std::pair<Point,unsigned>.
//the unsigned integer is incremented at each call to
//operator()
struct Auto_count : public CGAL::cpp98::unary_function<const Point&,std::pair<Point,unsigned> >{
  mutable unsigned i;
  Auto_count() : i(0){}
  std::pair<Point,unsigned> operator()(const Point& p) const {
    return std::make_pair(p,i++);
  }
};

int main()
{
  std::vector<Point> points;
  points.push_back(Point(0,0,0));
  points.push_back(Point(1,0,0));
  points.push_back(Point(0,1,0));
  points.push_back(Point(0,0,1));
  points.push_back(Point(2,2,2));
  points.push_back(Point(-1,0,1));

  Triangulation T( boost::make_transform_iterator(points.begin(),Auto_count()),
                   boost::make_transform_iterator(points.end(),  Auto_count() )  );

  assert( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Triangulation::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    if( points[ vit->info() ] != vit->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
