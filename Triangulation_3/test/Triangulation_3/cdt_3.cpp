#define CGAL_DEBUG_CDT_3 1
#define CGAL_TRIANGULATION_CHECK_EXPENSIVE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Conforming_Delaunay_triangulation_3.h>
#include <CGAL/draw_triangulation_3.h>

#include <vector>
#include <cassert>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Triangulation_data_structure_3<
  CGAL::Conforming_Delaunay_triangulation_vertex_base_3<K>,
  CGAL::Delaunay_triangulation_cell_base_3<K> >                 Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                  Delaunay;
typedef Delaunay::Point                                         Point;
using Point_3 = K::Point_3;

int main()
{
  std::cerr.precision(17);
  std::cout.precision(17);

  auto test1 = []() {
    const std::vector<Point_3> points = { { -2,  0,  0 },
                                          {  2,  0,  0 },
                                          {  0,  1, -1 },
                                          {  0, -1, -1 },
                                          {  0,  0,  1 },
                                          {  -10, -10, -10  },
                                          {  -10, 10, -10   },
                                          {  10, 10, -10    },
                                          {  10, -10, -10   },
                                          {  -10, -10, 10   },
                                          {  -10, 10, 10    },
                                          {  10, 10, 10     },
                                          {  10, -10, 10    },
    };
    std::vector<Delaunay::Vertex_handle> vertices;
    vertices.reserve(points.size());
    CGAL::Conforming_Delaunay_triangulation_3<Delaunay> cdt;
    for(auto p: points) vertices.push_back(cdt.insert(p));
    Delaunay::Cell_handle c;
    assert( cdt.is_valid() );
    assert(cdt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
    assert(cdt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

    Delaunay::Cell_handle ch;
    int li, lj;
    assert(!cdt.is_edge(vertices[0], vertices[1], ch, li, lj));
    cdt.insert_constrained_edge(vertices[0], vertices[1]);
    // cdt.insert_constrained_edge(vertices[5], vertices[1]);
    cdt.insert_constrained_edge(vertices[5], vertices[11]);
    return cdt.is_conforming() ? 0 : 1;
  };
  CGAL_USE(test1);

  auto test2 = []() {
    CGAL::Conforming_Delaunay_triangulation_3<Delaunay> cdt;

    std::ifstream input("clusters2.edg");
    if(!input) return 1;
    int n;
    input >> n;
    std::cerr << n << " lines in the file\n";
    while(n-- > 0) {
      double x, y;
      input >> x >> y;
      auto v1 = cdt.insert({x, y, 0});
        CGAL_assertion(cdt.is_conforming());
      auto v3 = cdt.insert({x, y, 1});
        CGAL_assertion(cdt.is_conforming());
      input >> x >> y;
      auto v2 = cdt.insert({x, y, 0});
        CGAL_assertion(cdt.is_conforming());
      auto v4 = cdt.insert({x, y, 1});
        CGAL_assertion(cdt.is_conforming());
      cdt.insert_constrained_edge(v1, v2);
        CGAL_assertion(cdt.is_conforming());
      cdt.insert_constrained_edge(v3, v4);
        CGAL_assertion(cdt.is_conforming());
    }
    cdt.insert({11, 28, 0});
      CGAL_assertion(cdt.is_conforming());
    cdt.insert({11, 33, 0});
      CGAL_assertion(cdt.is_conforming());
    // assert(cdt.is_edge(vertices[0], vertices[1], ch, li, lj));
    for(auto v: cdt.finite_vertex_handles()) {
      std::cout << "Point ( " << v->point() << " )\n";
      std::cout << "  on " << v->nb_of_incident_constraints << " constraint(s): "
                << v->c_id << "\n";
    }
    CGAL::draw(static_cast<Delaunay&>(cdt), "CDT_3", true);
    return cdt.is_conforming() ? 0 : 1;
  };

  return test1() + test2();
}
