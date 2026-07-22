#include <type_traits>
#include <iterator>
#include <CGAL/type_traits.h>

// Add C++20 utilities to C++14
template <class Iterator>
using iter_value_t = std::remove_reference_t<decltype(*std::declval<Iterator>())>;

template<class Range>
using range_value_t = iter_value_t<CGAL::cpp20::remove_cvref_t<decltype(std::begin(std::declval<Range>()))>>;

template< class T, class U >
constexpr bool is_same_v = std::is_same<T, U>::value;

template< class T, class U >
constexpr bool is_convertible_v = std::is_convertible<T, U>::value;


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

using Triangulation = CGAL::Triangulation_3<CGAL::Exact_predicates_inexact_constructions_kernel>;

using Vertex_handle =  Triangulation::Vertex_handle;

Triangulation t;

// Tds::vertex_handles() -> Vertex_handle
static_assert(is_same_v<range_value_t<decltype(t.tds().vertex_handles())>, const Vertex_handle>);

// Triangulation::all_vertex_handles() -> Vertex_handle
static_assert(is_same_v<range_value_t<decltype(t.all_vertex_handles())>, const Vertex_handle>);

// Triangulation::finite_vertex_handles() -> convertible to Vertex_handle
static_assert(is_convertible_v<decltype(*std::begin(t.finite_vertex_handles())), Vertex_handle>);

// But is it equal to Vertex_handle
static_assert(is_same_v<range_value_t<decltype(t.finite_vertex_handles())>, const Vertex_handle>);

int main()
{
  Vertex_handle v_inf = t.infinite_vertex();
  for(auto v: t.finite_vertex_handles()) {
    Vertex_handle v2 = v;
    CGAL_USE(v2);
    if(v == v_inf) return 1;
  }
  for (auto v : t.all_vertex_handles()) {
      Vertex_handle v2 = v;
      CGAL_USE(v2);
      if (v == v_inf) std::cout << "found inf" << std::endl;
  }
  return 0;
}
