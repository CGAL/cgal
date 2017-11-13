#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
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
void test_np(const NamedParameters& np)
{
// test values
  assert(   get_param(np,CGAL::internal_np::vertex_index).v == 0 );
  assert(   get_param(np,CGAL::internal_np::face_index).v == 2 );
  assert(   get_param(np,CGAL::internal_np::vertex_point).v == 3 );
  assert(   get_param(np,CGAL::internal_np::edge_is_constrained).v == 12 );
  assert(   get_param(np,CGAL::internal_np::density_control_factor).v == 14 );
  assert(   get_param(np,CGAL::internal_np::use_delaunay_triangulation).v == 15 );
  assert(   get_param(np,CGAL::internal_np::fairing_continuity).v == 16 );
  assert(   get_param(np,CGAL::internal_np::sparse_linear_solver).v == 17 );
  assert(   get_param(np,CGAL::internal_np::weight_calculator).v == 18 );

//test types
  check_same_type<0>( get_param(np,CGAL::internal_np::vertex_index) );
  check_same_type<2>( get_param(np,CGAL::internal_np::face_index) );
  check_same_type<3>( get_param(np,CGAL::internal_np::vertex_point) );
  check_same_type<12>( get_param(np,CGAL::internal_np::edge_is_constrained) );
  //
  check_same_type<14>( get_param(np,CGAL::internal_np::density_control_factor) );
  check_same_type<15>( get_param(np,CGAL::internal_np::use_delaunay_triangulation) );
  check_same_type<16>( get_param(np,CGAL::internal_np::fairing_continuity) );
  check_same_type<17>( get_param(np,CGAL::internal_np::sparse_linear_solver) );
  check_same_type<18>( get_param(np,CGAL::internal_np::weight_calculator) );
}

int main()
{
  test_np(
  // start with inherited named params
    CGAL::Polygon_mesh_processing::parameters::
                      vertex_index_map(A<0>(0)).
                      face_index_map(A<2>(2)).
                      vertex_point_map(A<3>(3)).
                      edge_is_constrained_map(A<12>(12)).
  // continue with PMP specific named params
                      density_control_factor(A<14>(14)).
                      use_delaunay_triangulation(A<15>(15)).
                      fairing_continuity(A<16>(16)).
                      sparse_linear_solver(A<17>(17)).
                      weight_calculator(A<18>(18))
  );
}
