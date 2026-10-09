#include <CGAL/Named_function_parameters.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <cassert>
#include <cstdlib>
#include <functional>
#include <type_traits>

namespace inp = CGAL::internal_np;
namespace params = CGAL::parameters;

template <int i>
using Static_int = std::integral_constant<int, i>;

void test_all_cgal_named_params() {
  struct A{};
  A a;
#define CGAL_add_named_parameter(X, Y, Z) \
  (void)params::Z(a).Z(a);
#include <CGAL/STL_Extension/internal/parameters_interface.h>
}

struct Non_copyable
{
  int value = 0;
  Non_copyable() = default;
  Non_copyable(const Non_copyable&) = delete;
};

template <int i, class T>
void check_same_type(T)
{
  static const bool b = std::is_same_v<Static_int<i>, T>;
  static_assert(b);
  assert(b);
}

void test_values_and_types()
{
  auto np = params::vertex_index_map(Static_int<0>())
                          .visitor(Static_int<1>());
  using params::get_parameter;
  using CGAL::parameter_or;

  // test values
  assert(get_parameter(np, inp::vertex_index).value == 0);
  assert(np.parameter(inp::vertex_index).value == 0);
  assert(get_parameter(np, inp::visitor).value == 1);
  assert(np.parameter(inp::visitor).value == 1);

  // test types
  check_same_type<0>(get_parameter(np, inp::vertex_index));
  check_same_type<0>(np.parameter(inp::vertex_index));
  check_same_type<0>(parameter_or(np, inp::vertex_index, Static_int<42>{}));
  check_same_type<0>(parameter_or<Static_int<42>>(np, inp::vertex_index));
  check_same_type<1>(get_parameter(np, inp::visitor));
  check_same_type<1>(np.parameter(inp::visitor));
  check_same_type<1>(parameter_or(np, inp::visitor, Static_int<42>{}));
  check_same_type<1>(parameter_or<Static_int<42>>(np, inp::visitor));

  auto v = parameter_or(np, inp::face_color_map, Static_int<42>{});
  assert(v.value == 42);
  assert(v.value == np.parameter_or(inp::face_color_map, Static_int<42>{}).value);

  auto v2 = parameter_or<Static_int<42>>(np, inp::face_color_map);
  assert(v2.value == 42);
  assert(v2.value == np.parameter_or<Static_int<42>>(inp::face_color_map).value);
}

void test_missing_parameters()
{
  const auto np = params::default_values();
  using NamedParameters = decltype(np);

  static_assert(!NamedParameters::has_parameter(inp::vertex_index));
  static_assert(params::is_default_parameter<NamedParameters, inp::vertex_index_t>::value);
  static_assert(std::is_same_v<decltype(np.parameter(inp::vertex_index)), inp::Param_not_found>);

  Static_int<4> fallback;
  Static_int<4>& result = CGAL::parameter_or(np, inp::vertex_index, fallback);
  assert(&result == &fallback);

  Static_int<4>& member_result = np.parameter_or(inp::vertex_index, fallback);
  assert(&member_result == &fallback);

  auto default_constructed = CGAL::parameter_or<Static_int<42>>(np, inp::vertex_index);
  assert(default_constructed.value == 42);
}

void test_compatibility_aliases()
{
  auto np = params::default_values()
    .seeds(Static_int<5>{})
    .time_limit(0.5)
    .max_iteration_number(9);

  static_assert(!params::is_default_parameter<decltype(np), inp::seeds_t>::value);
  static_assert(!params::is_default_parameter<decltype(np), inp::maximum_running_time_t>::value);
  static_assert(!params::is_default_parameter<decltype(np), inp::number_of_iterations_t>::value);
  assert(params::get_parameter(np, inp::seeds).value == 5);
  assert(params::get_parameter(np, inp::maximum_running_time) == 0.5);
  assert(params::get_parameter(np, inp::number_of_iterations) == 9);
}

void test_no_copyable()
{
  Non_copyable b;
  auto np = params::visitor(b);
  using NamedParameters = decltype(np);
  using NP_type = typename inp::Get_param<typename NamedParameters::base,inp::visitor_t>::type;
  static_assert(std::is_same_v<NP_type, std::reference_wrapper<const Non_copyable>>);

  const Static_int<4>& a = params::choose_parameter(
    params::get_parameter_reference(np, inp::edge_index), Static_int<4>());
  assert(a.value == 4);
}

void test_references()
{
  Non_copyable b;
  auto v = Static_int<0>();
  auto np = params::visitor(std::ref(b))
                                .vertex_point_map(b)
                                .vertex_index_map(v)
                                .face_index_map(std::cref(b));
  using NamedParameters = decltype(np);
  using Default_type = Static_int<2>;
  Default_type default_value;

  // std::reference_wrapper
  using Visitor_reference_type =
      typename inp::Lookup_named_param_def<inp::visitor_t, NamedParameters, Default_type>::reference;
  static_assert(std::is_same_v<Non_copyable&, Visitor_reference_type>);
  Visitor_reference_type vis_ref =
      params::choose_parameter(params::get_parameter_reference(np, inp::visitor), default_value);
  CGAL_USE(vis_ref);
  assert(&vis_ref == &b);
  vis_ref.value = 42;
  assert(b.value == 42);

  // std::reference_wrapper of const
  using FIM_reference_type =
      typename inp::Lookup_named_param_def<inp::face_index_t, NamedParameters, Default_type>::reference;
  static_assert(std::is_same_v<const Non_copyable&, FIM_reference_type>);
  FIM_reference_type fim_ref =
      params::choose_parameter(params::get_parameter_reference(np, inp::face_index), default_value);
  CGAL_USE(fim_ref);
  assert(&fim_ref == &b);

  // non-copyable
  using VPM_reference_type =
      typename inp::Lookup_named_param_def<inp::vertex_point_t, NamedParameters, Default_type>::reference;
  static_assert(std::is_same_v<const Non_copyable&, VPM_reference_type>);
  VPM_reference_type vpm_ref =
      params::choose_parameter(params::get_parameter_reference(np, inp::vertex_point), default_value);
  CGAL_USE(vpm_ref);
  assert(&vpm_ref == &b);

  // passed by copy
  using VIM_reference_type =
      typename inp::Lookup_named_param_def<inp::vertex_index_t, NamedParameters, Default_type>::reference;
  static_assert(std::is_same_v<Static_int<0>, VIM_reference_type>);
  VIM_reference_type vim_ref =
      params::choose_parameter(params::get_parameter_reference(np, inp::vertex_index), default_value);
  CGAL_USE(vim_ref);
  assert(&vim_ref != &v);

  // default
  using EIM_reference_type =
      typename inp::Lookup_named_param_def<inp::edge_index_t, NamedParameters, Default_type>::reference;
  static_assert(std::is_same_v<Default_type&, EIM_reference_type>);
  EIM_reference_type eim_ref =
      params::choose_parameter(params::get_parameter_reference(np, inp::edge_index), default_value);
  assert(&eim_ref == &default_value);
}


void test_ref_only_parameters()
{
  int i = 42;
  auto np_ref_only_i = params::weights(i).weights(i);
  auto& ref_i = np_ref_only_i.parameter_ref(CGAL::internal_np::weights_param_t{});
  assert(&ref_i == &i);

  const int ci = 43;
  auto np_cref_only_ci = params::image(ci).image(ci);
  auto& ref_ci = np_cref_only_ci.parameter_ref(CGAL::internal_np::image_3_param_t{});
  assert(&ref_ci == &ci);

  auto np_ref_only_ci = params::weights(ci).weights(ci);
  auto& ref_ci2 = np_ref_only_ci.parameter_ref(CGAL::internal_np::weights_param_t{});
  assert(&ref_ci2 == &ci);

  auto np_cref_only_2 = params::image(2*3).image(2*3);
  auto& ref_2 = np_cref_only_2.parameter_ref(CGAL::internal_np::image_3_param_t{});
  static_assert(std::is_const_v<std::remove_reference_t<decltype(ref_2)>>);
  assert(np_cref_only_2.parameter_ref(CGAL::internal_np::image_3_param_t{}) == 6);
}

void test_authorized_options()
{
  auto np_ok1 = CGAL::parameters::vertex_point_map(0).edge_index_map(2).face_index_map(3);
  auto np_ok2 = CGAL::parameters::face_index_map(3).vertex_point_map(0).edge_index_map(4);
  constexpr auto np_ko = CGAL::parameters::face_index_map(3).vertex_point_map(0).edge_index_map(4).halfedge_index_map(4);
  auto np_permissive = np_ko.do_not_check_allowed_np(true);
  auto np_default = CGAL::parameters::default_values();

  static_assert(decltype(np_ok1)::number_of_parameters == 3);
  static_assert(decltype(np_ok2)::number_of_parameters == 3);
  static_assert(decltype(np_ko)::number_of_parameters == 4);
  static_assert(decltype(np_permissive)::number_of_parameters == 5);
  static_assert(decltype(np_default)::number_of_parameters == 1);

  static_assert(CGAL::parameters::authorized_options<decltype(np_default),
                                                     CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>());
  static_assert(CGAL::parameters::authorized_options<decltype(np_ok1),
                                                     CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>());
  static_assert(CGAL::parameters::authorized_options<decltype(np_ok2),
                                                     CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>());
  static_assert(!CGAL::parameters::authorized_options<decltype(np_ko),
                                                      CGAL::internal_np::vertex_point_t,
                                                      CGAL::internal_np::edge_index_t,
                                                      CGAL::internal_np::face_index_t>());
  static_assert(CGAL::parameters::authorized_options<decltype(np_permissive),
                                                     CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>());
  static_assert(CGAL::parameters::authorized_options<decltype(np_ko.do_not_check_allowed_np(true)),
                                                     CGAL::internal_np::vertex_point_t,
                                                     CGAL::internal_np::edge_index_t,
                                                     CGAL::internal_np::face_index_t>());

  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS(np_ok1, vertex_point_t, edge_index_t, face_index_t);
  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS(np_ok2, vertex_point_t, edge_index_t, face_index_t);
}

int main()
{
  test_all_cgal_named_params();
  test_values_and_types();

  test_missing_parameters();

  test_compatibility_aliases();

  test_no_copyable();

  test_references();

  test_ref_only_parameters();

  // test that, in case of duplicates, the last parameter value is kept
  auto np = params::visitor(1).visitor(2);
  static_assert(decltype(np)::has_parameter(inp::visitor));
  assert(params::get_parameter(np, inp::visitor) == 2);
  assert(CGAL::parameter_or(np, inp::visitor, 3) == 2);

  auto d = CGAL::parameters::default_values();
  static_assert(std::is_same_v<decltype(d), CGAL::parameters::Default_named_parameters>);
#ifndef CGAL_NO_DEPRECATED_CODE
  auto d1 = CGAL::parameters::all_default();
  static_assert(std::is_same_v<decltype(d1), CGAL::parameters::Default_named_parameters>);
  auto d2 = CGAL::Polygon_mesh_processing::parameters::all_default();
  static_assert(std::is_same_v<decltype(d2), CGAL::parameters::Default_named_parameters>);
#endif

  test_authorized_options();

  return EXIT_SUCCESS;
}
