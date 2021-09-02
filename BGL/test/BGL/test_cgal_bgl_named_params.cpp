#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

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
  static const bool b = boost::is_same< A<i>, T >::value;
  CGAL_static_assertion(b);
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
void test_no_copyable(const NamedParameters&)
{
  typedef typename inp::Get_param<typename NamedParameters::base,inp::visitor_t>::type NP_type;
  static_assert( boost::is_same<NP_type,std::reference_wrapper<const B> > ::value );
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

  return EXIT_SUCCESS;
}
