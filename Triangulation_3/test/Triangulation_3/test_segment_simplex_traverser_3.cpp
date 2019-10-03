#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <array>

#include <CGAL/Random.h>

//#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >                DT;

typedef DT::Cell_handle                                         Cell_handle;

typedef CGAL::Triangulation_segment_simplex_iterator_3<DT>      Simplex_traverser;

template <typename Big_tuple>
bool test(const DT dt, const Big_tuple& tuple);

int main(int, char* [])
{
  const std::vector<Point_3> points = { { -2,  0,  0 },
                                        {  2,  0,  0 }, 
                                        {  0,  1,  -1 },
                                        {  0, -1,  -1 },
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
  std::vector<DT::Vertex_handle> vertices;
  vertices.reserve(points.size());
  DT dt;
  for(auto p: points) vertices.push_back(dt.insert(p));
  DT::Cell_handle c;
  assert( dt.is_valid() );
  assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

  std::cerr << dt.number_of_finite_cells() << '\n';

  const std::vector < std::tuple<Point_3, Point_3, std::array<unsigned, 4>>> queries = {
      {{-1, 0,  0}, { 2, 0,  0}, {1, 0, 1, 2}}, // CFCV
      // {{ 2, 0,  0}, {-1, 0,  0}, {1, 0, 1, 2}}, // reverse: assertion
      {{-2, 0,  0}, { 2, 0,  0}, {2, 0, 1, 2}}, // VCFCV
      {{ 2, 0,  0}, {-2, 0,  0}, {2, 0, 1, 2}}, // reverse case: VCFCV
      {{-3, 0,  0}, { 3, 0,  0}, {2, 0, 1, 4}}, // bug
      {{-2, 0,  0}, { 2, 2, -2}, {2, 1, 0, 1}}, // bug
      {{ 2, 2, -2}, {-2, 0,  0}, {2, 1, 0, 1}}, // bug
  };
  bool ok = true;
  for(const auto& tuple: queries) {
    if(!test(dt, tuple)) ok = false;
  }
  std::cout << "Done (" << queries.size() << " queries)\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Big_tuple>
bool test(const DT dt, const Big_tuple& tuple) {
  bool result = true;
  using std::get;
  const auto& p1 = get<0>(tuple);
  const auto& p2 = get<1>(tuple);
  const auto& expected_results = get<2>(tuple);

  std::cout << "\n#\n# Query segment: ( " << p1 << " , "
            << p2 << " )\n#\n";
  Simplex_traverser st(dt, p1, p2);
  
  unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;
  unsigned int nb_collinear = 0;

  // Count the number of finite cells traversed.
  unsigned int inf = 0, fin = 0;
  for (; st != st.end(); ++st) {
    if (Cell_handle(st) != Cell_handle() && dt.is_infinite(st))
      ++inf;
    else {
      ++fin;
      if (st.is_facet())
        ++nb_facets;
      else if (st.is_edge())
        ++nb_edges;
      else if (st.is_vertex())
        ++nb_vertex;
      else if (st.is_cell())
        ++nb_cells;

      if (st.is_collinear()) {
        ++nb_collinear;
        std::cout << "collinear\n";
      }

      switch (st->dimension()) {
      case 2: {
        std::cout << "facet " << std::endl;
        DT::Facet f = *st;
        std::cout << "  ( " << f.first->vertex((f.second + 1) & 3)->point()
                  << "  " << f.first->vertex((f.second + 2) & 3)->point()
                  << "  " << f.first->vertex((f.second + 3) & 3)->point()
                  << " )\n";
        break;
      }
      case 1: {
        std::cout << "edge " << std::endl;
        DT::Edge e = *st;
        std::cout << "  ( " << e.first->vertex(e.second)->point() << "  "
                  << e.first->vertex(e.third)->point() << " )\n";
        break;
      }
      case 0: {
        std::cout << "vertex " << std::endl;
        DT::Vertex_handle v = *st;
        std::cout << "  ( " << v->point() << " )\n";
        break;
      }
      case 3: {
        std::cout << "cell \n  ( ";
        DT::Cell_handle ch = *st;
        for (int i = 0; i < 4; ++i)
          std::cout << ch->vertex(i)->point() << "  ";
        std::cout << " )\n";
        break;
      }
      default:
        CGAL_assume(false);
      } // end switch
    }
      }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "While traversing from " << st.source()
                << " to " << st.target() << std::endl;
      std::cout << "\tinfinite cells : " << inf << std::endl;
      std::cout << "\tfinite cells   : " << fin << std::endl;
      std::cout << "\tfacets   : " << nb_facets << std::endl;
      std::cout << "\tedges    : " << nb_edges << std::endl;
      std::cout << "\tvertices : " << nb_vertex << std::endl;
      std::cout << std::endl << std::endl;
#endif

    auto check_expected = [](unsigned value, unsigned expected) {
      if(value != expected) {
        std::cout << "\t  ERROR: expected " << expected << '\n';
        return false;
      }
      return true;
    };
    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    result = result && check_expected(nb_cells, expected_results[3]);
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    result = result && check_expected(nb_facets, expected_results[2]);
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    result = result && check_expected(nb_edges, expected_results[1]);
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;
    result = result && check_expected(nb_vertex, expected_results[0]);
    std::cout << "\tnb collinear : " << nb_collinear << std::endl;
    return result;
}
