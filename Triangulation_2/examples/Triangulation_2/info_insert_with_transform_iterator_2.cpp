#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Vertex_handle                                     Vertex_handle;

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
  points.push_back(Point(0,0));
  points.push_back(Point(1,0));
  points.push_back(Point(0,1));
  points.push_back(Point(4,10));
  points.push_back(Point(2,2));
  points.push_back(Point(-1,0));


  Delaunay T;
  T.insert( boost::make_transform_iterator(points.begin(),Auto_count()),
            boost::make_transform_iterator(points.end(),  Auto_count() )  );

  CGAL_assertion( T.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (Vertex_handle v : T.finite_vertex_handles())
    if( points[ v->info() ] != v->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "OK" << std::endl;

  return 0;
}
