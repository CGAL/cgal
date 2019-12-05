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
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location>  Delaunay;
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
    Delaunay dt;
    for(auto p: points) vertices.push_back(dt.insert(p));
    Delaunay::Cell_handle c;
    assert( dt.is_valid() );
    assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
    assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

    std::cerr << dt.number_of_vertices() << '\n';
    std::cerr << dt.number_of_finite_cells() << '\n';
    CGAL::Triangulation_conformer_3<Delaunay> conformer(dt);
    Delaunay::Cell_handle ch;
    int li, lj;
    assert(!dt.is_edge(vertices[0], vertices[1], ch, li, lj));
    conformer.insert_constrained_edge(vertices[0], vertices[1]);
    // conformer.insert_constrained_edge(vertices[5], vertices[1]);
    conformer.restore_Delaunay();
    conformer.insert_constrained_edge(vertices[5], vertices[11]);
    conformer.restore_Delaunay();
    return 0;
  };
  CGAL_USE(test1);

  auto test2 = []() {
    Delaunay dt;
    CGAL::Triangulation_conformer_3<Delaunay> conformer(dt);

    std::ifstream input("clusters2.edg");
    if(!input) return 1;
    int n;
    input >> n;
    std::cerr << n << " lines in the file\n";
    while(n-- > 0) {
      double x, y;
      input >> x >> y;
      auto v1 = dt.insert({x, y, 0});
      auto v3 = dt.insert({x, y, 1});
      input >> x >> y;
      auto v2 = dt.insert({x, y, 0});
      auto v4 = dt.insert({x, y, 1});
      conformer.insert_constrained_edge(v1, v2);
      conformer.insert_constrained_edge(v3, v4);
    }
    conformer.restore_Delaunay();
    std::cerr << dt.number_of_vertices() << '\n';
    std::cerr << dt.number_of_finite_cells() << '\n';
    // assert(dt.is_edge(vertices[0], vertices[1], ch, li, lj));
    CGAL::draw(dt, "CDT_3", true);
    return 0;
  };

  return test1() + test2();
}
