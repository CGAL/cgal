

#include <fstream>
#include <cassert>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Simple_cartesian<double> K1;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
typedef CGAL::Triangulation_3<K1>      Tr1;
typedef CGAL::Triangulation_3<K2>      Tr2;


template <typename T1, typename T2>
struct Update_vertex
{
  typedef typename T1::Vertex                  V1;
  typedef typename T2::Vertex                  V2;
  typedef typename T2::Point                   Point;

  V2 operator()(const V1&)
  {
    return V2();
  }

  void operator()(const V1& v1, V2& v2)
  {
    CGAL::Cartesian_converter<K1, K2> c;
    v2.set_point(Point(c(v1.point())));
  }
}; // end struct Update_vertex

struct Update_cell {
  template <typename C1, typename C2>
  void operator()(const C1&, C2&) {}
}; // end struct Update_cell
int main()
{
  // construction from a list of points :
  std::list<Tr1::Point> L;
  L.push_back(Tr1::Point(0,0,0));
  L.push_back(Tr1::Point(1,0,0));
  L.push_back(Tr1::Point(0,1,0));
  L.push_back(Tr1::Point(0,0,1));
  Tr1 T1(L.begin(), L.end());
  std::ofstream out("tr");
  out << T1;
  out.close();

  Tr2 T2;
  std::ifstream in("tr");
  T2.file_input<Tr1,Update_vertex<Tr1, Tr2>, Update_cell>(in);
  in.close();
  assert(T2.is_valid());
  Tr2::Point_iterator pit = T2.points_begin();
  assert(*(pit)++ == Tr2::Point(0,0,0));
  assert(*(pit)++ == Tr2::Point(1,0,0));
  assert(*(pit)++ == Tr2::Point(0,1,0));
  assert(*(pit)++ == Tr2::Point(0,0,1));


  std::cout << "done" << std::endl;
  return 0;
}
