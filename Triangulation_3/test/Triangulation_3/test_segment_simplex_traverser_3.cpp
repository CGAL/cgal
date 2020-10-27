#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <array>

#include <CGAL/Random.h>

#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3<Kernel> DT;
typedef DT::Cell_handle                        Cell_handle;
typedef DT::Segment_simplex_iterator           Segment_simplex_iterator;

bool test(const DT& dt,
          const std::pair<Point_3, Point_3>& query,
          const std::array<unsigned, 4>& expected_result);

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
  Cell_handle c;
  assert(dt.is_valid());
  assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

  std::cerr << dt.number_of_finite_cells() << '\n';

  const std::vector < std::pair<Point_3, Point_3>> queries = {
      {{-1, 0,  0}, { 1, 0,  0}}, // CFC
      {{-1, 0,  0}, { 2, 0,  0}}, // CFCV
      {{ 2, 0,  0}, {-1, 0,  0}}, // reverse
      {{-2, 0,  0}, { 2, 0,  0}}, // VCFCV
      {{ 2, 0,  0}, {-2, 0,  0}}, // reverse case: VCFCV
      {{-3, 0,  0}, { 3, 0,  0}}, // FVCFCVF
      {{-2, 0,  0}, { 2, 2, -2}}, // VEVF
      {{ 2, 2, -2}, {-2, 0,  0}} // reverse case: FVEV
  };

  const std::vector< std::array<unsigned, 4> > expected_results ={
      {0, 0, 1, 2}, // CFC
      {1, 0, 1, 2}, // CFCV
      {1, 0, 1, 2}, // reverse
      {2, 0, 1, 2}, // VCFCV
      {2, 0, 1, 2}, // reverse case: VCFCV
      {2, 0, 3, 2}, // FVCFCVF
      {2, 1, 1, 0}, // VEVF
      {2, 1, 1, 0} // reverse case: FVEV
  };

  bool ok = true;
  for(std::size_t i=0; i<queries.size(); ++i) {
    if(!test(dt, queries[i], expected_results[i])) ok = false;
  }
  std::cout << "Done (" << queries.size() << " queries)\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

bool test(const DT& dt,
          const std::pair<Point_3, Point_3>& query,
          const std::array<unsigned, 4>& expected_result)
{
  bool result = true;
  using std::get;
  const auto& p1 = query.first;
  const auto& p2 = query.second;
  unsigned int nb_cells = 0, nb_collinear = 0;
  unsigned int nb_facets = 0, nb_edges = 0, nb_vertex = 0;

  std::cout << "\n#\n# Query segment: ( " << p1 << " , "
            << p2 << " )\n#\n";
  Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(p1, p2);
  Segment_simplex_iterator stend = dt.segment_traverser_simplices_end();

  // Count the number of finite cells traversed.
  unsigned int inf = 0, fin = 0;
  for (; st != stend; ++st)
  {
    if (st->dimension() == 3
      && dt.is_infinite(Cell_handle(*st)))
      ++inf;
    else {
      ++fin;

      switch (st->dimension()) {
      case 2: {
        ++nb_facets;
        std::cout << "facet " << std::endl;
        DT::Facet f = *st;
        std::cout << "  ( " << f.first->vertex((f.second + 1) & 3)->point()
                  << "  " << f.first->vertex((f.second + 2) & 3)->point()
                  << "  " << f.first->vertex((f.second + 3) & 3)->point()
                  << " )\n";
        break;
      }
      case 1: {
        ++nb_edges;
        std::cout << "edge " << std::endl;
        DT::Edge e = *st;
        std::cout << "  ( " << e.first->vertex(e.second)->point() << "  "
                  << e.first->vertex(e.third)->point() << " )\n";
        break;
      }
      case 0: {
        ++nb_vertex;
        std::cout << "vertex " << std::endl;
        DT::Vertex_handle v = *st;
        std::cout << "  ( " << v->point() << " )\n";
        break;
      }
      case 3: {
        ++nb_cells;
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
    result = result && check_expected(nb_cells, expected_result[3]);
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    result = result && check_expected(nb_facets, expected_result[2]);
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    result = result && check_expected(nb_edges, expected_result[1]);
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;
    result = result && check_expected(nb_vertex, expected_result[0]);
    std::cout << "\tnb collinear : " << nb_collinear << std::endl;

    return result;
}
