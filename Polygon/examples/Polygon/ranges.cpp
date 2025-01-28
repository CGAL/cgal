#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef K::Point_2                                          Point_2;
typedef K::Segment_2                                        Segment_2;


int main()
{
  // create a polygon and put some points in it
  Polygon_2 p;

  p.push_back(Point_2(0,0));
  p.push_back(Point_2(4,0));
  p.push_back(Point_2(4,4));


  for(const Point_2& p : p.vertices()){
    std::cout << p << std::endl;
  }

  // As the range is not light weight, we have to use a reference
  const Polygon_2::Vertices& range = p.vertices();

  for(auto it = range.begin(); it!= range.end(); ++it){
    std::cout << *it << std::endl;
  }

  for(const Segment_2& e  : p.edges()){
    std::cout << e << std::endl;
  }

  return EXIT_SUCCESS;
}
