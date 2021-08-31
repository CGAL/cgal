#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

#include <cstdlib>

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
  using CGAL::parameters::get_parameter;

  // test values
  assert(get_parameter(np, CGAL::internal_np::vertex_index).v == 0);
  assert(get_parameter(np, CGAL::internal_np::visitor).v == 1);

  // test types
  check_same_type<0>(get_parameter(np, CGAL::internal_np::vertex_index));
  check_same_type<1>(get_parameter(np, CGAL::internal_np::visitor));
}

template<class NamedParameters>
void test_no_copyable(const NamedParameters&)
{
  typedef typename CGAL::internal_np::Get_param<typename NamedParameters::base,CGAL::internal_np::visitor_t>::type NP_type;
  static_assert( boost::is_same<NP_type,std::reference_wrapper<const B> > ::value );
}

int main()
{
  test_values_and_types(CGAL::parameters::vertex_index_map(A<0>(0)).visitor(A<1>(1)));

  B b;
  test_no_copyable(CGAL::parameters::visitor(b));
  return EXIT_SUCCESS;
}
