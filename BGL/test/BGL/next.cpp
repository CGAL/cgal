
#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;
typedef boost::graph_traits<SM>::halfedge_descriptor halfedge_descriptor;


namespace Toto {

  struct Graph {};

  struct Descriptor{};

  Descriptor next(const Descriptor& d, const Graph&)
  {
    return d;
  }

}

namespace std {

  template <>
  struct iterator_traits<Toto::Descriptor> {

  };
};


namespace Nested {
  void
  gnu()
  {
    SM sm;
    CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,0),sm);

    halfedge_descriptor hd = *(halfedges(sm).first);
    hd = next(hd,sm);
    std::cout << hd;
  }

  void
  gnats()
  {
    Toto::Graph g;
    Toto::Descriptor d;
    d = next(d,g);
  }

}

int main()
{
  Nested::gnu();

  Nested::gnats();

  return 0;
}
