#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

template <int i>
struct A{
  A(int v):v(v){}
  int v;
};

template <int i, class T>
void check_same_type(T)
{
  static const bool b = boost::is_same< A<i>, T >::value;
  CGAL_static_assertion(b);
  assert(b);
}

template<class NamedParameters>
void test(const NamedParameters& np)
{
// test values
  assert(   get_param(np,boost::vertex_index).v == 0 );
  assert(   get_param(np,boost::halfedge_index).v == 1 );
  assert(   get_param(np,boost::face_index).v == 2 );
  assert(   get_param(np,boost::vertex_point).v == 3 );
  assert(   get_param(np,CGAL::vertex_is_fixed).v == 4 );
  assert(   get_param(np,boost::edge_index).v == 5 );
  assert(   get_param(np,boost::graph_visitor).v == 6 );
  assert(   get_param(np,CGAL::set_cache_policy).v == 7 );
  assert(   get_param(np,CGAL::get_cost_policy).v == 8 );
  assert(   get_param(np,CGAL::get_cost_policy_params).v == 9 );
  assert(   get_param(np,CGAL::get_placement_policy).v == 10 );
  assert(   get_param(np,CGAL::get_placement_policy_params).v == 11 );
  assert(   get_param(np,CGAL::edge_is_constrained).v == 12 );
  assert(   get_param(np,CGAL::edge_is_constrained_params).v == 13 );

//test types
  check_same_type<0>( get_param(np,boost::vertex_index) );
  check_same_type<1>( get_param(np,boost::halfedge_index) );
  check_same_type<2>( get_param(np,boost::face_index) );
  check_same_type<3>( get_param(np,boost::vertex_point) );
  check_same_type<4>( get_param(np,CGAL::vertex_is_fixed) );
  check_same_type<5>( get_param(np,boost::edge_index) );
  check_same_type<6>( get_param(np,boost::graph_visitor) );
  check_same_type<7>( get_param(np,CGAL::set_cache_policy) );
  check_same_type<8>( get_param(np,CGAL::get_cost_policy) );
  check_same_type<9>( get_param(np,CGAL::get_cost_policy_params) );
  check_same_type<10>( get_param(np,CGAL::get_placement_policy) );
  check_same_type<11>( get_param(np,CGAL::get_placement_policy_params) );
  check_same_type<12>( get_param(np,CGAL::edge_is_constrained) );
  check_same_type<13>( get_param(np,CGAL::edge_is_constrained_params) );
}

int main()
{
  test(
    CGAL::parameters::vertex_index_map(A<0>(0)).
                      halfedge_index_map(A<1>(1)).
                      face_index_map(A<2>(2)).
                      vertex_point_map(A<3>(3)).
                      vertex_is_fixed_map(A<4>(4)).
                      edge_index_map(A<5>(5)).
                      visitor(A<6>(6)).
                      set_cache(A<7>(7)).
                      get_cost(A<8>(8)).
                      get_cost_params(A<9>(9)).
                      get_placement(A<10>(10)).
                      get_placement_params(A<11>(11)).
                      edge_is_constrained_map(A<12>(12)).
                      edge_is_constrained_map_params(A<13>(13))
  );
}
