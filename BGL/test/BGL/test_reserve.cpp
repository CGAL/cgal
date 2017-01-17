#include "test_Prefix.h"

#include <CGAL/boost/graph/helpers.h>

template<typename Mesh>
void test() {
  Mesh m;
  CGAL::reserve(m, 6, 12, 7);
}

struct MyGraph{};
struct MyGraphWithReserve{
  bool OK;
  MyGraphWithReserve() 
    : OK(false) 
  {}
  void reserve(int,int,int)
  {
    OK=true;
  }
};

namespace boost
{
  template<>
  struct graph_traits<MyGraph>
  {
    typedef int vertices_size_type;
    typedef int edges_size_type;
    typedef int faces_size_type;
  };
  template<>
  struct graph_traits<MyGraphWithReserve>
    : graph_traits<MyGraph>
  {};
}

int main()
{
  test<SM>();
  test<Polyhedron>();
#if defined(CGAL_USE_OPENMESH)
  test<OMesh>();
#endif
  test<MyGraph>();
  MyGraphWithReserve g;
  assert(!g.OK);
  CGAL::reserve(g,1,2,3);
  assert(g.OK);
  return 0;
}
