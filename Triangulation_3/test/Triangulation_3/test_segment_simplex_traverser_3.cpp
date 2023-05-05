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

// Define the structure.
typedef CGAL::Base_with_time_stamp<CGAL::Triangulation_vertex_base_3<Kernel>> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> DT;
typedef DT::Cell_handle                        Cell_handle;
typedef DT::Edge                               Edge;
typedef DT::Facet                              Facet;
typedef DT::Vertex_handle                      Vertex_handle;
typedef DT::Segment_simplex_iterator           Segment_simplex_iterator;

// a function to insert without spatial sorting
template <typename Point_it>
void insert(DT& dt, Point_it first, Point_it end) {
  for(; first != end; ++first) {
    dt.insert(*first);
  }
}
// trick to avoid a conflict with CGAL-5.6
//  -> specialization of Output_rep for `CC_iterator` *and* `My`
struct My {};
template <class DSC, bool Const >
class CGAL::Output_rep<CGAL::internal::CC_iterator<DSC, Const>, My > {
protected:
  using CC_iterator = CGAL::internal::CC_iterator<DSC, Const>;
  using Compact_container = typename CC_iterator::CC;
  using Time_stamper = typename Compact_container::Time_stamper;
  CC_iterator it;
public:
  Output_rep( const CC_iterator it) : it(it) {}
  std::ostream& operator()( std::ostream& out) const {
    return out << '#' << it->time_stamp();
  }
};
auto display_vert(Vertex_handle v) {
  std::stringstream os;
  os.precision(17);
  os << CGAL::IO::oformat(v, My()) << "=(" << v->point() << ")";
  return os.str();
};
template <typename Simplex>
struct Debug_simplex {
  Simplex simplex;

  template<typename CharT, typename Traits>
  friend
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const Debug_simplex& d) {
    auto&& simplex = d.simplex;
    switch(simplex.dimension()) {
      case 0: {
        os << "   vertex " << display_vert(static_cast<Vertex_handle>(simplex));
        break;
      }
      case 1: {
        const auto [c, index1, index2] = static_cast<Edge>(simplex);
        os << "   egde "
           << display_vert(c->vertex(index1)) << " - "
           << display_vert(c->vertex(index2));
        break;
      }
      case 2: {
        const auto [c, index] = static_cast<Facet>(simplex);
        os << "   facet "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 0))) << " - "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 1))) << " - "
           << display_vert(c->vertex(DT::vertex_triple_index(index, 2)));
        break;
      }
      case 3: {
        const auto c = static_cast<Cell_handle>(simplex);
        os << "   cell "
           << display_vert(c->vertex(0)) << " - "
           << display_vert(c->vertex(1)) << " - "
           << display_vert(c->vertex(2)) << " - "
           << display_vert(c->vertex(3));
        break;
      }
      default: CGAL_assume(false);
    }
    return os;
  };
};
template <typename Simplex>
auto debug_simplex(Simplex simplex) {
  return Debug_simplex<Simplex>{simplex};
}

static const std::vector<Point_3> bbox_points =
{
  {  -10.1, -10, -10  },
  {  -10.2, 10, -10   },
  {  10.3, 10, -10    },
  {  10.4, -10, -10   },
  {  -10.5, -10, 10   },
  {  -10.6, 10, 10    },
  {  10.7, 10, 10     },
  {  10.8, -10, 10    },
  };

DT dt;
std::string result_string;

auto visit_simplex = [](auto s) {
  auto d = s.dimension();
  if(3 == d && dt.is_infinite(static_cast<Cell_handle>(s))) {
    result_string += 'I';
  } else {
    result_string += std::to_string(d);
  }
  std::cout << debug_simplex(s) << '\n';
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

  for(auto s: dt.segment_traverser_simplices(v[0], v[3])) {
    visit_simplex(s);
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

bool test_a_simple_tetrahedron() {
  std::cerr << "## test_a_simple_tetrahedron()\n";
  bool ok = true;
  auto test = [&](Point_3 a, Point_3 b, bool with_bbox, std::string expected_result) {
    bool exception_thrown = false;
    dt.clear();
    dt.insert({0, 0, 0});
    dt.insert({1, 0, 0});
    dt.insert({0, 1, 0});
    dt.insert({0, 0, 1});
    if(with_bbox) insert(dt, bbox_points.begin(), bbox_points.end());
    result_string.clear();
    std::cerr << "### Case " << expected_result;
    if(with_bbox) std::cerr << " with bbox";
    std::cerr << '\n';
    try {
      for(auto s: dt.segment_traverser_simplices(a, b)) {
        visit_simplex(s);
      }
    } catch(const CGAL::Assertion_exception& e) {
      CGAL::get_static_warning_handler()("Assertion", e.expression().c_str(),
                                         e.filename().c_str(),
                                         e.line_number(),
                                         e.message().c_str());
      exception_thrown = true;
    }
    if(result_string != expected_result) {
      std::cerr << "test_a_simple_tetrahedron failed\n";
      std::cerr << "  result_string is " << result_string << " instead of "
                << expected_result << '\n';
      ok = false;
    }
    if(exception_thrown) {
      std::cerr << "test_a_simple_tetrahedron failed\n";
      std::cerr << "  exception thrown\n";
      ok = false;
    }
  };

  test({ 0,  0,  0}, { 1,  0,  0}, false, "010");
  test({ 0,  0,  0}, { 2,  0,  0}, false, "010I");
  test({-1,  0,  0}, { 2,  0,  0}, false, "I010I");
  test({ 0,  0,  0}, {.5, .5,  0}, false, "021");
  test({ 0,  0,  0}, { 1,  1,  0}, false, "021I");
  test({-1, -1,  0}, { 1,  1,  0}, false, "I021I");

  test({ 0,  0,  0}, { 1,  0,  0}, true, "010");
  test({ 0,  0,  0}, { 2,  0,  0}, true, "0103");
  test({-1,  0,  0}, { 2,  0,  0}, true, "30103");
  test({ 0,  0,  0}, {.5, .5,  0}, true, "021");
  test({ 0,  0,  0}, { 1,  1,  0}, true, "0213");
  test({-1, -1,  0}, { 1,  1,  0}, true, "30213");

  return ok;
}

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
  dt.clear();
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
  ok = test_vfefv() && ok;
  ok = test_vfefv(true) && ok;
  ok = test_a_simple_tetrahedron() && ok;
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

      visit_simplex(*st);

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
