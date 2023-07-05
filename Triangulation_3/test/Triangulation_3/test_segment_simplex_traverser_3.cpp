//#define CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/IO/io.h>

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
typedef Kernel::Segment_3                                       Segment_3;
typedef Kernel::Triangle_3                                      Triangle_3;
typedef Kernel::Tetrahedron_3                                   Tetrahedron_3;

// Define the structure.
typedef CGAL::Base_with_time_stamp<CGAL::Triangulation_vertex_base_3<Kernel>> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> DT;
typedef DT::Cell_handle                        Cell_handle;
typedef DT::Edge                               Edge;
typedef DT::Facet                              Facet;
typedef DT::Vertex_handle                      Vertex_handle;
typedef DT::Simplex                            Simplex;
typedef DT::Segment_simplex_iterator           Segment_simplex_iterator;

#include "test_triangulation_simplex_3_debug.h"

// a function to insert without spatial sorting
template <typename Point_it>
void insert(DT& dt, Point_it first, Point_it end) {
  for(; first != end; ++first) {
    dt.insert(*first);
  }
}

static const std::vector<Point_3> bbox_points =
{
  { -10.1, -10, -10.08  },
  { -10.2,  10, -10.07  },
  {  10.3,  10, -10.06  },
  {  10.4, -10, -10.05  },
  { -10.5, -10,  10.04  },
  { -10.6,  10,  10.03  },
  {  10.7,  10,  10.02  },
  {  10.8, -10,  10.01  },
  };

DT dt;
std::string result_string;

bool reverse_sort_vertex_handles(Vertex_handle v1, Vertex_handle v2) {
  return v1.operator->() > v2.operator->();
};

auto vertices_of_simplex(Simplex simplex) -> std::array<Vertex_handle, 4> {
  std::array<Vertex_handle, 4> vertices = { Vertex_handle{}, Vertex_handle{}, Vertex_handle{}, Vertex_handle{} };
  switch(simplex.dimension()) {
    case 0: {
      vertices[0] = static_cast<Vertex_handle>(simplex);
      break;
    }
    case 1: {
      const auto [c, index1, index2] = static_cast<Edge>(simplex);
      vertices[0] = c->vertex(index1);
      vertices[1] = c->vertex(index2);
      break;
    }
    case 2: {
      const auto [c, index] = static_cast<Facet>(simplex);
      vertices[0] = c->vertex(DT::vertex_triple_index(index, 0));
      vertices[1] = c->vertex(DT::vertex_triple_index(index, 1));
      vertices[2] = c->vertex(DT::vertex_triple_index(index, 2));
      break;
    }
    case 3: {
      const auto c = static_cast<Cell_handle>(simplex);
      vertices[0] = c->vertex(0);
      vertices[1] = c->vertex(1);
      vertices[2] = c->vertex(2);
      vertices[3] = c->vertex(3);
      break;
    }
    default: CGAL_unreachable();
  }
  std::sort(vertices.begin(), vertices.end(), reverse_sort_vertex_handles);
  for(int i = 0; i < 4; ++i) {
    assert((i <= simplex.dimension()) == (vertices[i] != Vertex_handle{}));
  }
  return vertices;
}

std::variant<Point_3, Segment_3, Triangle_3, Tetrahedron_3> get_simplex_geometry(Simplex simplex) {
  switch(simplex.dimension()) {
    case 0: {
      return static_cast<Vertex_handle>(simplex)->point();
    }
    case 1: {
      const auto [c, index1, index2] = static_cast<Edge>(simplex);
      return Segment_3(c->vertex(index1)->point(), c->vertex(index2)->point());
    }
    case 2: {
      const auto [c, index] = static_cast<Facet>(simplex);
      return Triangle_3(c->vertex(DT::vertex_triple_index(index, 0))->point(),
                        c->vertex(DT::vertex_triple_index(index, 1))->point(),
                        c->vertex(DT::vertex_triple_index(index, 2))->point());
    }
    case 3: {
      const auto c = static_cast<Cell_handle>(simplex);
      return Tetrahedron_3(c->vertex(0)->point(),
                           c->vertex(1)->point(),
                           c->vertex(2)->point(),
                           c->vertex(3)->point());
    }
    default: CGAL_unreachable();
  }
}

void visit_simplex(Point_3 a, Point_3 b, Simplex s, std::optional<Simplex> previous_simplex_optional) {
  auto d = s.dimension();
  if(3 == d && dt.is_infinite(static_cast<Cell_handle>(s))) {
    result_string += 'I';
  } else {
    result_string += std::to_string(d);
  }
  std::clog << debug_simplex(s) << '\n';
  const bool does_intersect_ab = (3 == d && dt.is_infinite(static_cast<Cell_handle>(s))) ||
    std::visit(
      [&](auto geometry) { return CGAL::do_intersect(Segment_3(a, b), geometry); },
      get_simplex_geometry(s));
  if(!does_intersect_ab) {
    CGAL_error_msg("the simplex does not intersect the query segment");
  }
  if(previous_simplex_optional) {
    // this block checks that consecutive simplices are incident
    using Set = std::array<Vertex_handle, 4>;
    Set prev_vertices = vertices_of_simplex(*previous_simplex_optional);
    Set s_vertices = vertices_of_simplex(s);
    if(previous_simplex_optional->dimension() < s.dimension()) {
      std::swap(prev_vertices, s_vertices);
      std::swap(*previous_simplex_optional, s);
    }
    if(!std::includes(
            prev_vertices.begin(), prev_vertices.begin() + 1 + previous_simplex_optional->dimension(),
            s_vertices.begin(),    s_vertices.begin() + 1 + s.dimension(),
            reverse_sort_vertex_handles))
    {
      CGAL_error_msg("consecutive simplices are not incident");
    }
  }
};

bool test_vfefv(bool with_bbox = false)
{
  std::cerr << "## test_vfefv(" << std::boolalpha << with_bbox << ")\n";
  result_string.clear();
  static const std::vector<Point_3> points =
  {
    { -1,  0,  0 },
    {  0,  1,  0 },
    {  0, -1,  0 },
    {  5,  0,  0 },
    {  6,  2,  2 },
    {  6, -2, -2 },
  };

  dt.clear();
  insert(dt, points.begin(), points.end());
  if(with_bbox) insert(dt, bbox_points.begin(), bbox_points.end());

  const auto v = dt.finite_vertex_handles().to<std::vector>();

  Cell_handle c; int i, j, k;
  assert(dt.is_facet(v[0], v[1], v[2], c, i, j, k));
  assert(dt.is_facet(v[1], v[2], v[3], c, i, j, k));
  assert(dt.is_cell (v[1], v[2], v[3], v[4], c));
  assert(dt.is_cell (v[1], v[2], v[3], v[5], c));

  std::optional<Simplex> previous{};
  for(auto s: dt.segment_traverser_simplices(v[0], v[3])) {
    visit_simplex(points[0], points[3], s, previous);
    previous = s;
  }
  static const std::string expected_result_string = "02120";
  bool ok = (result_string == expected_result_string);
  if(!ok) {
    std::cerr << "test_vfefv failed\n";
    std::cerr << "  result_string is " << result_string << " instead of "
              << expected_result_string << '\n';
  }
  return ok;
}

bool test_a_simple_tetrahedron(const std::vector<Point_3>& points) {

  DT dt2;
  for(const auto& p: points) dt2.insert(p);

  bool ok = true;
  auto test = [&](Point_3 a, Point_3 b, std::string expected_result) {
    // This test function calls `do_test` with four configurations:
    //  - with [ab] and [ba],
    //  - and with or without a bbox around the central tetrahedron.
    dt = dt2;
    auto do_with_or_without_bbox = [&](Point_3 a, Point_3 b, bool with_bbox, std::string expected_result) {
      auto do_it = [&](auto from, auto to) {
        bool exception_thrown = false;
        result_string.clear();
        try {
          std::optional<Simplex> previous_simplex;
          for(auto s: dt.segment_traverser_simplices(from, to)) {
            visit_simplex(a, b, s, previous_simplex);
            previous_simplex = s;
          }
        } catch(const CGAL::Assertion_exception& e) {
          CGAL::get_static_warning_handler()("Assertion", e.expression().c_str(),
                                             e.filename().c_str(),
                                             e.line_number(),
                                             e.message().c_str());
          exception_thrown = true;
        }
        if(result_string != expected_result || exception_thrown) {
          std::cerr << "test_a_simple_tetrahedron failed on case " << expected_result
                    << (with_bbox ? " with bbox\n" : "\n");
          ok = false;
        }
        if(result_string != expected_result) {
          std::cerr << "  result_string is " << result_string << " instead of "
                    << expected_result << '\n';
        }
        if(exception_thrown) {
          std::cerr << "  exception thrown\n";
        }
      }; // end do_it

      std::clog << "### Case " << expected_result;
      if(with_bbox) std::clog << " with bbox";
      std::clog << '\n';
      std::clog << "from (" << a << ") to (" << b << ")\n";
      do_it(a, b);

      // then re-test using vertex handles, if possible
      Vertex_handle va{};
      Vertex_handle vb{};
      DT::Locate_type lt;
      int i, j;
      auto c = dt.locate(a, lt, i, j);
      if(lt == DT::VERTEX) va = c->vertex(i);
      c = dt.locate(b, lt, i, j);
      if(lt == DT::VERTEX) vb = c->vertex(i);
      if(va != Vertex_handle{} && vb != Vertex_handle{}) {
        std::clog << "from vertex" << display_vert(va) << " to vertex" << display_vert(vb) << ")\n";
        do_it(va, vb);
      };
    }; // end do_with_or_without_bbox
    std::string expected_result_reversed{expected_result.rbegin(), expected_result.rend()};
    do_with_or_without_bbox(a, b, false, expected_result);
    do_with_or_without_bbox(b, a, false, expected_result_reversed);
    std::replace(expected_result.begin(), expected_result.end(), 'I', '3');
    std::replace(expected_result_reversed.begin(), expected_result_reversed.end(), 'I', '3');
    insert(dt, bbox_points.begin(), bbox_points.end());
    do_with_or_without_bbox(a, b, true, expected_result);
    do_with_or_without_bbox(b, a, true, expected_result_reversed);
  }; // end test() lambda

  // [010] queries entering by a vertex and exiting by the other vertex of an incident edge,
  // on the line (x,0,0)
  test({ 0,  0,  0}, {.5,  0,  0},  "01");
  test({ 0,  0,  0}, { 1,  0,  0},  "010");
  test({ 0,  0,  0}, { 2,  0,  0},  "010I");
  test({-1,  0,  0}, { 2,  0,  0}, "I010I");
  test({-1,  0,  0}, { 1,  0,  0}, "I010");
  test({-1,  0,  0}, {.5,  0,  0}, "I01");

  // x [020] is not possible, because that would have passed through the edge -> see [010]

  // [021] queries entering by a vertex and exiting by the opposite edge of a same face,
  // on the line (x,x,0) (y==x)
  test({ 0,  0,  0}, {.2, .2,  0},  "02");
  test({ 0,  0,  0}, {.5, .5,  0},  "021");
  test({ 0,  0,  0}, { 1,  1,  0},  "021I");
  test({-1, -1,  0}, { 1,  1,  0}, "I021I");
  test({-1, -1,  0}, {.5, .5,  0}, "I021");
  test({-1, -1,  0}, {.2, .2,  0}, "I02");

  // x [030] (entering by a vertex and exiting by a non-adjacent vertex) is not possible
  //   because that would have passed by the common edge -> see [010]

  // x [031] (entering by a vertex and exiting by a non-incident edge), is not possible
  //   because that would have passed by the face -> see [021]

  // [032] queries entering by a vertex and exiting by the opposite facet,
  // on the line x==y==0.25-0.25z
  test({  0,   0,   1}, { .25,    .25,     .25 },  "03");
  test({  0,   0,   1}, { .25,    .25,     0   },  "032");
  test({  0,   0,   1}, { .5,     .5,     -.1  },  "032I");
  test({-.25,-.25,  2}, { .28125, .28125, -.125}, "I032I");
  test({-.25,-.25,  2}, { .25,    .25,     0   }, "I032");
  test({-.25,-.25,  2}, { .125,   .125,    .5  }, "I03");

  // [121] queries entering by an edge and exiting by an edge of the same face,
  // on the line (x,.5,0)
  test({0,   .5,  0}, {.2,  .5,  0},  "12");
  test({0,   .5,  0}, {.5 , .5,  0},  "121");
  test({0,   .5,  0}, {.6 , .5,  0},  "121I");
  test({-.1, .5,  0}, {.6 , .5,  0}, "I121I");
  test({-.1, .5,  0}, {.5 , .5,  0}, "I121");
  test({-.1, .5,  0}, {.25, .5,  0}, "I12");

  // x [130] (entering by an edge and exiting by a non-incident vertex) is not possible
  //   because that would have passed by the common face -> see [120] ([021] in reverse)

  // [131] entering by an edge and exiting by the opposite edge in the cell,
  // on the line x==y==0.5-z, also known as (.5-z, .5-z, z)
  test({  0,     0,     .5  }, { .25,  .25,   .25 },   "13");
  test({  0,     0,     .5  }, { .5,   .5,    0   },   "131");
  test({  0,     0,     .5  }, { .625, .625, -.125},   "131I");
  test({ -.125, -.125,  .625}, { .625, .625, -.125},  "I131I");
  test({ -.125, -.125,  .625}, { .5 ,  .5 ,   0   },  "I131");
  test({ -.125, -.125,  .625}, { .25,  .25,   .25 },  "I13");

  // [132] queries entering by an edge and exiting by a facet, on the line (x, .25-x, x)
  test({  0, .25,  0}, { .20, .05,  .20},  "13");
  test({  0, .25,  0}, { .25,   0,  .25},  "132");
  test({  0, .25,  0}, { .5 ,-.25,  .5 },  "132I");
  test({-.5, .75,-.5}, { .5 ,-.25,  .5 }, "I132I");
  test({-.5, .75,-.5}, { .25,   0,  .25}, "I132");
  test({-.5, .75,-.5}, { .20, .05,  .20}, "I13");

  // [232] queries entering by a facet and exiting by a facet, on the line (x,.5-x,.2)
  test({ 0,   .5,  .2}, {.2,   .3,  .2},  "23");
  test({ 0,   .5,  .2}, {.5,  0,    .2},  "232");
  test({ 0,   .5,  .2}, {.55, -.05, .2},  "232I");
  test({-.05, .45, .2}, {.55, -.05, .2}, "I232I");
  test({-.05, .45, .2}, {.5,  0,    .2}, "I232");
  test({-.05, .45, .2}, {.2,   .3,  .2}, "I23");

  // special case: queries stay in a single simplex
  test({ -.125, -.125,  .625}, { -.125, -.125,  .6251},    "I");
  test({ .25, .25,  .25},      { .20, .25,  .25},      "3");
  test({ 0,   .5,  .2},        { 0,   .4,  .2},        "2");
  test({ 0,   .5,  .0},        { 0,   .6,  .0},        "1");

  return ok;
}

bool test_a_simple_tetrahedron() {
  bool ok = true;
  std::cout << "## test_a_simple_tetrahedron()\n"
            << R"#(
This test uses with a trivial tetrahedron, and is launched with all the
24 permutations of the four vertices. There are 7 test cases, and for each of
them 6 different segments are tested:
- The segment starts at the incoming boundary of the tetrahedron and ends
  inside.
- The segment starts at the incoming boundary of the tetrahedron and ends
  at the outgoing boundary.
- The segment starts at the incoming boundary of the tetrahedron, goes out
  and beyond.
- The segment starts before the tetrahedron, goes through it, and comes out.
- The segment starts before the tetrahedron and ends at the outgoing boundary.
- The segment starts before the tetrahedron and ends inside.

For each of those query segments, 4 tests are performed:
  - with/without 8 extra vertices in the triangulation, forming a bounding
    box,
  - and in the direct and reverse direction.

In total, 4032 tests are performed...


)#";
  std::cout.flush();
  std::vector<Point_3> points = {
      {0., 0., 0.},
      {1., 0., 0.},
      {0., 1., 0.},
      {0., 0., 1.},
  };
  std::sort(points.begin(), points.end());
  do {
    std::cout << "### new permutation of the four points\n";
    for(const auto& p: points) std::cout << "    " << p << '\n';
    std::cout << std::endl;
    ok = test_a_simple_tetrahedron(points) && ok;
  }
  while (std::next_permutation(points.begin(), points.end()));
  return ok;
}

bool test(const DT& dt,
          const std::pair<Point_3, Point_3>& query,
          const std::array<unsigned, 4>& expected_result);


int main(int, char* [])
{
  std::cerr.precision(17);
  std::clog.precision(17);
  std::cout.precision(17);

  bool ok = true;
  ok = test_a_simple_tetrahedron() && ok;

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
  dt.clear();
  for(auto p: points) vertices.push_back(dt.insert(p));
  Cell_handle c;
  assert(dt.is_valid());
  assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

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

  for(std::size_t i=0; i<queries.size(); ++i) {
    ok = test(dt, queries[i], expected_results[i]) && ok;
  }
  std::cout << "Done (" << queries.size() << " queries)\n";

  ok = test_vfefv() && ok;
  ok = test_vfefv(true) && ok;
  if(ok) {
    std::cout << "All tests passed\n";
    return EXIT_SUCCESS;
  } else {
    std::cout << "Some tests failed\n";
    return EXIT_FAILURE;
  }
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
  std::optional<Simplex> previous;
  for (; st != stend; ++st)
  {
    if (st->dimension() == 3
      && dt.is_infinite(Cell_handle(*st)))
      ++inf;
    else {
      ++fin;

      visit_simplex(p1, p2, *st, previous);
      previous = *st;

      switch (st->dimension()) {
      case 2: ++nb_facets; break;
      case 1: ++nb_edges; break;
      case 0: ++nb_vertex; break;
      case 3: ++nb_cells; break;
      default: CGAL_unreachable();
      }
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
