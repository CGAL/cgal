#include <CGAL/Named_function_parameters.h>
#include <CGAL/assertions.h>
#include <type_traits>

#include <cstdlib>

namespace inp = CGAL::internal_np;
namespace params = CGAL::parameters;

template <int i>
struct A
{
  A(int v):v(v){}
  int v;
};

struct B
{
  B(){}
  B(const B&) = delete;
};

template <int i, class T>
void check_same_type(T)
{
  static const bool b = std::is_same< A<i>, T >::value;
  static_assert(b);
  assert(b);
}

template<class NamedParameters>
void test_values_and_types(const NamedParameters& np)
{
  using params::get_parameter;

  // test values
  assert(get_parameter(np, inp::vertex_index).v == 0);
  assert(get_parameter(np, inp::visitor).v == 1);

  // test types
  check_same_type<0>(get_parameter(np, inp::vertex_index));
  check_same_type<1>(get_parameter(np, inp::visitor));
}

template<class NamedParameters>
void test_no_copyable(const NamedParameters& np)
{
  typedef typename inp::Get_param<typename NamedParameters::base,inp::visitor_t>::type NP_type;
  static_assert(std::is_same<NP_type,std::reference_wrapper<const B> >::value);

  const A<4>& a  = params::choose_parameter(params::get_parameter_reference(np, inp::edge_index), A<4>(4));
  assert(a.v==4);
}

template <class NamedParameters>
void test_references(const NamedParameters& np)
{
  typedef A<2> Default_type;
  Default_type default_value(2);

  // std::reference_wrapper
  typedef typename inp::Lookup_named_param_def<inp::visitor_t, NamedParameters, Default_type>::reference Visitor_reference_type;
  static_assert(std::is_same<B&, Visitor_reference_type>::value);
  Visitor_reference_type vis_ref = params::choose_parameter(params::get_parameter_reference(np, inp::visitor), default_value);
  CGAL_USE(vis_ref);

  // std::reference_wrapper of const
  typedef typename inp::Lookup_named_param_def<inp::face_index_t, NamedParameters, Default_type>::reference FIM_reference_type;
  static_assert(std::is_same<const B&, FIM_reference_type>::value);
  FIM_reference_type fim_ref = params::choose_parameter(params::get_parameter_reference(np, inp::face_index), default_value);
  CGAL_USE(fim_ref);

  // non-copyable
  typedef typename inp::Lookup_named_param_def<inp::vertex_point_t, NamedParameters, Default_type>::reference VPM_reference_type;
  static_assert(std::is_same<const B&, VPM_reference_type>::value);
  VPM_reference_type vpm_ref = params::choose_parameter(params::get_parameter_reference(np, inp::vertex_point), default_value);
  CGAL_USE(vpm_ref);

  // passed by copy
  typedef typename inp::Lookup_named_param_def<inp::vertex_index_t, NamedParameters, Default_type>::reference VIM_reference_type;
  static_assert(std::is_same<A<0>, VIM_reference_type>::value);
  VIM_reference_type vim_ref = params::choose_parameter(params::get_parameter_reference(np, inp::vertex_index), default_value);
  CGAL_USE(vim_ref);

  // default
  typedef typename inp::Lookup_named_param_def<inp::edge_index_t, NamedParameters, Default_type>::reference EIM_reference_type;
  static_assert(std::is_same<Default_type&, EIM_reference_type>::value);
  EIM_reference_type eim_ref = params::choose_parameter(params::get_parameter_reference(np, inp::edge_index), default_value);
  assert(&eim_ref==&default_value);
}

void test_authorized_options()
{
  auto np_ok1 = CGAL::parameters::vertex_point_map(0).edge_index_map(2).face_index_map(3);
  auto np_ok2 = CGAL::parameters::face_index_map(3).vertex_point_map(0).edge_index_map(4);
  auto np_ko = CGAL::parameters::face_index_map(3).vertex_point_map(0).edge_index_map(4).halfedge_index_map(4);
  auto np_default = CGAL::parameters::default_values();


  static_assert(CGAL::parameters::authorized_options<CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>(np_default));
  static_assert(CGAL::parameters::authorized_options<CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>(np_ok1));
  static_assert(CGAL::parameters::authorized_options<CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>(np_ok2));
  static_assert(!CGAL::parameters::authorized_options<CGAL::internal_np::vertex_point_t,
                                                      CGAL::internal_np::edge_index_t,
                                                      CGAL::internal_np::face_index_t>(np_ko));

  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS(np_ok1, vertex_point_t, edge_index_t, face_index_t);
  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS(np_ok2, vertex_point_t, edge_index_t, face_index_t);
}

int main()
{
  test_values_and_types(params::vertex_index_map(A<0>(0)).visitor(A<1>(1)));

  B b;
  test_no_copyable(params::visitor(b));

  test_references(params::visitor(std::ref(b))
                         .vertex_point_map(b)
                         .vertex_index_map(A<0>(0))
                         .face_index_map(std::reference_wrapper<const B>(b))
  );

  auto d = CGAL::parameters::default_values();
  static_assert(std::is_same<decltype(d),CGAL::parameters::Default_named_parameters>::value);
#ifndef CGAL_NO_DEPRECATED_CODE
  auto d1 = CGAL::parameters::all_default();
  static_assert(std::is_same<decltype(d1),CGAL::parameters::Default_named_parameters>::value);
  auto d2 = CGAL::Polygon_mesh_processing::parameters::all_default();
  static_assert(std::is_same<decltype(d2),CGAL::parameters::Default_named_parameters>::value);
#endif

  test_authorized_options();

  return EXIT_SUCCESS;
}
