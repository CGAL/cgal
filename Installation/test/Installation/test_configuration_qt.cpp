#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/IO/io.h>
#include <fstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>                            Triangulation;
typedef Triangulation::Point                                Point;

template <typename A, typename B>
typename CGAL::Coercion_traits<A,B>::Type
binary_func(const A& a , const B& b){
  typedef CGAL::Coercion_traits<A,B> CT;
  static_assert(CT::Are_explicit_interoperable::value);
  typename CT::Cast cast;
  return cast(a)*cast(b);
}

int main(int argc, char**) {
  std::cout<< binary_func(double(3), int(5)) << std::endl;
  std::cout<< binary_func(int(3), double(5)) << std::endl;
  std::ifstream in("data/triangulation_prog1.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;
  Triangulation t;
  t.insert(begin, end);
  if(argc == 3) // do not test Qt6 at runtime
    CGAL::draw(t);
  std::cout<<"OK."<<std::endl;
  return EXIT_SUCCESS;
}
