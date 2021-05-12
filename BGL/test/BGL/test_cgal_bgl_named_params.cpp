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
  using CGAL::parameters::get_parameter;

  // Test values

    // Named parameters that we use in CGAL
  assert(get_parameter(np, CGAL::internal_np::vertex_index).v == 0);
  assert(get_parameter(np, CGAL::internal_np::visitor).v == 1);
  assert(get_parameter(np, CGAL::internal_np::vertex_point).v == 2);
  assert(get_parameter(np, CGAL::internal_np::halfedge_index).v == 3);
  assert(get_parameter(np, CGAL::internal_np::edge_index).v == 4);
  assert(get_parameter(np, CGAL::internal_np::face_index).v == 5);

  assert(get_parameter(np, CGAL::internal_np::edge_is_constrained).v == 6);
  assert(get_parameter(np, CGAL::internal_np::first_index).v == 7);
  assert(get_parameter(np, CGAL::internal_np::number_of_iterations).v == 8);

  assert(get_parameter(np, CGAL::internal_np::METIS_options).v == 800000001);
  assert(get_parameter(np, CGAL::internal_np::vertex_partition_id).v == 800000002);
  assert(get_parameter(np, CGAL::internal_np::face_partition_id).v == 800000003);

  assert(get_parameter(np, CGAL::internal_np::vertex_to_vertex_output_iterator).v == 800000004);
  assert(get_parameter(np, CGAL::internal_np::halfedge_to_halfedge_output_iterator).v == 800000005);
  assert(get_parameter(np, CGAL::internal_np::face_to_face_output_iterator).v == 800000006);

  assert(get_parameter(np, CGAL::internal_np::vertex_to_vertex_map).v == 800000007);
  assert(get_parameter(np, CGAL::internal_np::halfedge_to_halfedge_map).v == 800000008);
  assert(get_parameter(np, CGAL::internal_np::face_to_face_map).v == 800000009);

  assert(get_parameter(np, CGAL::internal_np::implementation_tag).v == 800000010);
  assert(get_parameter(np, CGAL::internal_np::prevent_unselection).v == 800000011);

  assert(get_parameter(np, CGAL::internal_np::stream_precision).v == 800000012);

    // Named parameters that we use in the package 'Mesh_3'
  assert(get_parameter(np, CGAL::internal_np::vertex_feature_degree).v == 9);

    // Named parameters used in the package 'Polygon Mesh Processing'
  assert(get_parameter(np, CGAL::internal_np::geom_traits).v == 10);
  assert(get_parameter(np, CGAL::internal_np::vertex_incident_patches).v == 11);
  assert(get_parameter(np, CGAL::internal_np::density_control_factor).v == 12);
  assert(get_parameter(np, CGAL::internal_np::use_delaunay_triangulation).v == 13);
  assert(get_parameter(np, CGAL::internal_np::use_2d_constrained_delaunay_triangulation).v == 4573);
  assert(get_parameter(np, CGAL::internal_np::fairing_continuity).v == 14);
  assert(get_parameter(np, CGAL::internal_np::sparse_linear_solver).v == 15);
  assert(get_parameter(np, CGAL::internal_np::number_of_relaxation_steps).v == 16);
  assert(get_parameter(np, CGAL::internal_np::protect_constraints).v == 17);
  assert(get_parameter(np, CGAL::internal_np::relax_constraints).v == 18);
  assert(get_parameter(np, CGAL::internal_np::collapse_constraints).v == 43);
  assert(get_parameter(np, CGAL::internal_np::vertex_is_constrained).v == 19);
  assert(get_parameter(np, CGAL::internal_np::face_patch).v == 20);
  assert(get_parameter(np, CGAL::internal_np::random_uniform_sampling).v == 21);
  assert(get_parameter(np, CGAL::internal_np::grid_sampling).v == 22);
  assert(get_parameter(np, CGAL::internal_np::monte_carlo_sampling).v == 23);
  assert(get_parameter(np, CGAL::internal_np::do_sample_edges).v == 24);
  assert(get_parameter(np, CGAL::internal_np::do_sample_vertices).v == 25);
  assert(get_parameter(np, CGAL::internal_np::do_sample_faces).v == 26);
  assert(get_parameter(np, CGAL::internal_np::number_of_points_on_faces).v == 27);
  assert(get_parameter(np, CGAL::internal_np::number_of_points_per_face).v == 28);
  assert(get_parameter(np, CGAL::internal_np::grid_spacing).v == 29);
  assert(get_parameter(np, CGAL::internal_np::number_of_points_per_edge).v == 30);
  assert(get_parameter(np, CGAL::internal_np::number_of_points_on_edges).v == 31);
  assert(get_parameter(np, CGAL::internal_np::nb_points_per_area_unit).v == 32);
  assert(get_parameter(np, CGAL::internal_np::nb_points_per_distance_unit).v == 33);
  assert(get_parameter(np, CGAL::internal_np::throw_on_self_intersection).v == 43);
  assert(get_parameter(np, CGAL::internal_np::clip_volume).v == 44);
  assert(get_parameter(np, CGAL::internal_np::use_compact_clipper).v == 45);
  assert(get_parameter(np, CGAL::internal_np::erase_all_duplicates).v == 48);
  assert(get_parameter(np, CGAL::internal_np::require_same_orientation).v == 49);
  assert(get_parameter(np, CGAL::internal_np::use_bool_op_to_clip_surface).v == 50);
  assert(get_parameter(np, CGAL::internal_np::face_size_map).v == 52);
  assert(get_parameter(np, CGAL::internal_np::use_angle_smoothing).v == 53);
  assert(get_parameter(np, CGAL::internal_np::use_area_smoothing).v == 54);
  assert(get_parameter(np, CGAL::internal_np::use_Delaunay_flips).v == 55);
  assert(get_parameter(np, CGAL::internal_np::use_safety_constraints).v == 56);
  assert(get_parameter(np, CGAL::internal_np::area_threshold).v == 57);
  assert(get_parameter(np, CGAL::internal_np::volume_threshold).v == 58);
  assert(get_parameter(np, CGAL::internal_np::snapping_tolerance).v == 59);
  assert(get_parameter(np, CGAL::internal_np::dry_run).v == 60);
  assert(get_parameter(np, CGAL::internal_np::do_lock_mesh).v == 61);
  assert(get_parameter(np, CGAL::internal_np::halfedges_keeper).v == 62);
  assert(get_parameter(np, CGAL::internal_np::do_simplify_border).v == 64);
  assert(get_parameter(np, CGAL::internal_np::do_not_modify).v == 65);
  assert(get_parameter(np, CGAL::internal_np::allow_self_intersections).v == 66);
  assert(get_parameter(np, CGAL::internal_np::polyhedral_envelope_epsilon).v == 67);
  assert(get_parameter(np, CGAL::internal_np::maximum_number_of_faces).v == 78910);
  assert(get_parameter(np, CGAL::internal_np::non_manifold_feature_map).v == 60);
  assert(get_parameter(np, CGAL::internal_np::filter).v == 61);
  assert(get_parameter(np, CGAL::internal_np::face_epsilon_map).v == 62);

    // Named parameters that we use in the package 'Surface Mesh Simplification'
  assert(get_parameter(np, CGAL::internal_np::get_cost_policy).v == 34);
  assert(get_parameter(np, CGAL::internal_np::get_placement_policy).v == 35);

  // Named parameters that we use in the package 'Optimal_bounding_box'
  assert(get_parameter(np, CGAL::internal_np::use_convex_hull).v == 63);

    // To-be-documented named parameters
  assert(get_parameter(np, CGAL::internal_np::face_normal).v == 36);
  assert(get_parameter(np, CGAL::internal_np::random_seed).v == 37);
  assert(get_parameter(np, CGAL::internal_np::do_project).v == 38);

    // Internal named parameters
  assert(get_parameter(np, CGAL::internal_np::weight_calculator).v == 39);
  assert(get_parameter(np, CGAL::internal_np::preserve_genus).v == 40);
  assert(get_parameter(np, CGAL::internal_np::verbosity_level).v == 41);
  assert(get_parameter(np, CGAL::internal_np::use_binary_mode).v == 51);
  assert(get_parameter(np, CGAL::internal_np::projection_functor).v == 42);
  assert(get_parameter(np, CGAL::internal_np::apply_per_connected_component).v == 46);
  assert(get_parameter(np, CGAL::internal_np::output_iterator).v == 47);

  // Test types

    // Named parameters that we use in CGAL
  check_same_type<0>(get_parameter(np, CGAL::internal_np::vertex_index));
  check_same_type<1>(get_parameter(np, CGAL::internal_np::visitor));
  check_same_type<2>(get_parameter(np, CGAL::internal_np::vertex_point));
  check_same_type<3>(get_parameter(np, CGAL::internal_np::halfedge_index));
  check_same_type<4>(get_parameter(np, CGAL::internal_np::edge_index));
  check_same_type<5>(get_parameter(np, CGAL::internal_np::face_index));

  check_same_type<6>(get_parameter(np, CGAL::internal_np::edge_is_constrained));
  check_same_type<7>(get_parameter(np, CGAL::internal_np::first_index));
  check_same_type<8>(get_parameter(np, CGAL::internal_np::number_of_iterations));

  check_same_type<800000001>(get_parameter(np, CGAL::internal_np::METIS_options));
  check_same_type<800000002>(get_parameter(np, CGAL::internal_np::vertex_partition_id));
  check_same_type<800000003>(get_parameter(np, CGAL::internal_np::face_partition_id));
  check_same_type<800000004>(get_parameter(np, CGAL::internal_np::vertex_to_vertex_output_iterator));
  check_same_type<800000005>(get_parameter(np, CGAL::internal_np::halfedge_to_halfedge_output_iterator));
  check_same_type<800000006>(get_parameter(np, CGAL::internal_np::face_to_face_output_iterator));
  check_same_type<800000007>(get_parameter(np, CGAL::internal_np::vertex_to_vertex_map));
  check_same_type<800000008>(get_parameter(np, CGAL::internal_np::halfedge_to_halfedge_map));
  check_same_type<800000009>(get_parameter(np, CGAL::internal_np::face_to_face_map));
  check_same_type<800000010>(get_parameter(np, CGAL::internal_np::implementation_tag));
  check_same_type<800000011>(get_parameter(np, CGAL::internal_np::prevent_unselection));
  check_same_type<800000012>(get_parameter(np, CGAL::internal_np::stream_precision));

    // Named parameters that we use in the package 'Mesh_3'
  check_same_type<9>(get_parameter(np, CGAL::internal_np::vertex_feature_degree));

    // Named parameters used in the package 'Polygon Mesh Processing'
  check_same_type<10>(get_parameter(np, CGAL::internal_np::geom_traits));
  check_same_type<11>(get_parameter(np, CGAL::internal_np::vertex_incident_patches));
  check_same_type<12>(get_parameter(np, CGAL::internal_np::density_control_factor));
  check_same_type<13>(get_parameter(np, CGAL::internal_np::use_delaunay_triangulation));
  check_same_type<4573>(get_parameter(np, CGAL::internal_np::use_2d_constrained_delaunay_triangulation));
  check_same_type<14>(get_parameter(np, CGAL::internal_np::fairing_continuity));
  check_same_type<15>(get_parameter(np, CGAL::internal_np::sparse_linear_solver));
  check_same_type<16>(get_parameter(np, CGAL::internal_np::number_of_relaxation_steps));
  check_same_type<17>(get_parameter(np, CGAL::internal_np::protect_constraints));
  check_same_type<18>(get_parameter(np, CGAL::internal_np::relax_constraints));
  check_same_type<43>(get_parameter(np, CGAL::internal_np::collapse_constraints));
  check_same_type<19>(get_parameter(np, CGAL::internal_np::vertex_is_constrained));
  check_same_type<20>(get_parameter(np, CGAL::internal_np::face_patch));
  check_same_type<21>(get_parameter(np, CGAL::internal_np::random_uniform_sampling));
  check_same_type<22>(get_parameter(np, CGAL::internal_np::grid_sampling));
  check_same_type<23>(get_parameter(np, CGAL::internal_np::monte_carlo_sampling));
  check_same_type<24>(get_parameter(np, CGAL::internal_np::do_sample_edges));
  check_same_type<25>(get_parameter(np, CGAL::internal_np::do_sample_vertices));
  check_same_type<26>(get_parameter(np, CGAL::internal_np::do_sample_faces));
  check_same_type<27>(get_parameter(np, CGAL::internal_np::number_of_points_on_faces));
  check_same_type<28>(get_parameter(np, CGAL::internal_np::number_of_points_per_face));
  check_same_type<29>(get_parameter(np, CGAL::internal_np::grid_spacing));
  check_same_type<30>(get_parameter(np, CGAL::internal_np::number_of_points_per_edge));
  check_same_type<31>(get_parameter(np, CGAL::internal_np::number_of_points_on_edges));
  check_same_type<32>(get_parameter(np, CGAL::internal_np::nb_points_per_area_unit));
  check_same_type<33>(get_parameter(np, CGAL::internal_np::nb_points_per_distance_unit));
  check_same_type<43>(get_parameter(np, CGAL::internal_np::throw_on_self_intersection));
  check_same_type<44>(get_parameter(np, CGAL::internal_np::clip_volume));
  check_same_type<45>(get_parameter(np, CGAL::internal_np::use_compact_clipper));
  check_same_type<48>(get_parameter(np, CGAL::internal_np::erase_all_duplicates));
  check_same_type<49>(get_parameter(np, CGAL::internal_np::require_same_orientation));
  check_same_type<50>(get_parameter(np, CGAL::internal_np::use_bool_op_to_clip_surface));
  check_same_type<52>(get_parameter(np, CGAL::internal_np::face_size_map));
  check_same_type<53>(get_parameter(np, CGAL::internal_np::use_angle_smoothing));
  check_same_type<54>(get_parameter(np, CGAL::internal_np::use_area_smoothing));
  check_same_type<55>(get_parameter(np, CGAL::internal_np::use_Delaunay_flips));
  check_same_type<56>(get_parameter(np, CGAL::internal_np::use_safety_constraints));
  check_same_type<65>(get_parameter(np, CGAL::internal_np::do_not_modify));
  check_same_type<66>(get_parameter(np, CGAL::internal_np::allow_self_intersections));
  check_same_type<67>(get_parameter(np, CGAL::internal_np::polyhedral_envelope_epsilon));

  check_same_type<12340>(get_parameter(np, CGAL::internal_np::do_self_intersection_tests));
  check_same_type<12341>(get_parameter(np, CGAL::internal_np::do_orientation_tests));
  check_same_type<12342>(get_parameter(np, CGAL::internal_np::error_codes));
  check_same_type<12343>(get_parameter(np, CGAL::internal_np::volume_inclusions));
  check_same_type<12344>(get_parameter(np, CGAL::internal_np::face_connected_component_map));
  check_same_type<12345>(get_parameter(np, CGAL::internal_np::connected_component_id_to_volume_id));
  check_same_type<12346>(get_parameter(np, CGAL::internal_np::is_cc_outward_oriented));
  check_same_type<12347>(get_parameter(np, CGAL::internal_np::intersecting_volume_pairs_output_iterator));
  check_same_type<12348>(get_parameter(np, CGAL::internal_np::i_used_as_a_predicate));
  check_same_type<12349>(get_parameter(np, CGAL::internal_np::nesting_levels));
  check_same_type<12350>(get_parameter(np, CGAL::internal_np::i_used_for_volume_orientation));

  check_same_type<57>(get_parameter(np, CGAL::internal_np::area_threshold));
  check_same_type<58>(get_parameter(np, CGAL::internal_np::volume_threshold));
  check_same_type<59>(get_parameter(np, CGAL::internal_np::snapping_tolerance));
  check_same_type<60>(get_parameter(np, CGAL::internal_np::dry_run));
  check_same_type<61>(get_parameter(np, CGAL::internal_np::do_lock_mesh));
  check_same_type<62>(get_parameter(np, CGAL::internal_np::halfedges_keeper));
  check_same_type<64>(get_parameter(np, CGAL::internal_np::do_simplify_border));
  check_same_type<78910>(get_parameter(np, CGAL::internal_np::maximum_number_of_faces));
  check_same_type<60>(get_parameter(np, CGAL::internal_np::non_manifold_feature_map));
  check_same_type<61>(get_parameter(np, CGAL::internal_np::filter));
  check_same_type<62>(get_parameter(np, CGAL::internal_np::face_epsilon_map));

    // Named parameters that we use in the package 'Surface Mesh Simplification'
  check_same_type<34>(get_parameter(np, CGAL::internal_np::get_cost_policy));
  check_same_type<35>(get_parameter(np, CGAL::internal_np::get_placement_policy));

  // Named parameters that we use in the package 'Optimal_bounding_box'
  check_same_type<63>(get_parameter(np, CGAL::internal_np::use_convex_hull));

    // To-be-documented named parameters
  check_same_type<36>(get_parameter(np, CGAL::internal_np::face_normal));
  check_same_type<37>(get_parameter(np, CGAL::internal_np::random_seed));
  check_same_type<38>(get_parameter(np, CGAL::internal_np::do_project));
  check_same_type<456>(get_parameter(np, CGAL::internal_np::algorithm));

    // Internal named parameters
  check_same_type<39>(get_parameter(np, CGAL::internal_np::weight_calculator));
  check_same_type<40>(get_parameter(np, CGAL::internal_np::preserve_genus));
  check_same_type<41>(get_parameter(np, CGAL::internal_np::verbosity_level));
  check_same_type<51>(get_parameter(np, CGAL::internal_np::use_binary_mode));
  check_same_type<42>(get_parameter(np, CGAL::internal_np::projection_functor));
  check_same_type<46>(get_parameter(np, CGAL::internal_np::apply_per_connected_component));
  check_same_type<47>(get_parameter(np, CGAL::internal_np::output_iterator));

  // Named parameters used in the package 'Point Set Processing'
  check_same_type<9000>(get_parameter(np, CGAL::internal_np::point_map));
  check_same_type<9001>(get_parameter(np, CGAL::internal_np::query_point_map));
  check_same_type<9002>(get_parameter(np, CGAL::internal_np::normal_map));
  check_same_type<9003>(get_parameter(np, CGAL::internal_np::diagonalize_traits));
  check_same_type<9004>(get_parameter(np, CGAL::internal_np::svd_traits));
  check_same_type<9005>(get_parameter(np, CGAL::internal_np::callback));
  check_same_type<9006>(get_parameter(np, CGAL::internal_np::sharpness_angle));
  check_same_type<9007>(get_parameter(np, CGAL::internal_np::edge_sensitivity));
  check_same_type<9008>(get_parameter(np, CGAL::internal_np::neighbor_radius));
  check_same_type<9009>(get_parameter(np, CGAL::internal_np::number_of_output_points));
  check_same_type<9010>(get_parameter(np, CGAL::internal_np::size));
  check_same_type<9011>(get_parameter(np, CGAL::internal_np::maximum_variation));
  check_same_type<9012>(get_parameter(np, CGAL::internal_np::degree_fitting));
  check_same_type<9013>(get_parameter(np, CGAL::internal_np::degree_monge));
  check_same_type<9014>(get_parameter(np, CGAL::internal_np::threshold_percent));
  check_same_type<9015>(get_parameter(np, CGAL::internal_np::threshold_distance));
  check_same_type<9016>(get_parameter(np, CGAL::internal_np::attraction_factor));
  check_same_type<9017>(get_parameter(np, CGAL::internal_np::plane_map));
  check_same_type<9018>(get_parameter(np, CGAL::internal_np::plane_index_map));
  check_same_type<9019>(get_parameter(np, CGAL::internal_np::select_percentage));
  check_same_type<9020>(get_parameter(np, CGAL::internal_np::require_uniform_sampling));
  check_same_type<9021>(get_parameter(np, CGAL::internal_np::point_is_constrained));
  check_same_type<9022>(get_parameter(np, CGAL::internal_np::number_of_samples));
  check_same_type<9023>(get_parameter(np, CGAL::internal_np::accuracy));
  check_same_type<9024>(get_parameter(np, CGAL::internal_np::maximum_running_time));
  check_same_type<9025>(get_parameter(np, CGAL::internal_np::overlap));
  check_same_type<9026>(get_parameter(np, CGAL::internal_np::transformation));
  check_same_type<9027>(get_parameter(np, CGAL::internal_np::point_set_filters));
  check_same_type<9028>(get_parameter(np, CGAL::internal_np::matcher));
  check_same_type<9029>(get_parameter(np, CGAL::internal_np::outlier_filters));
  check_same_type<9030>(get_parameter(np, CGAL::internal_np::error_minimizer));
  check_same_type<9031>(get_parameter(np, CGAL::internal_np::transformation_checkers));
  check_same_type<9032>(get_parameter(np, CGAL::internal_np::inspector));
  check_same_type<9033>(get_parameter(np, CGAL::internal_np::logger));
  check_same_type<9034>(get_parameter(np, CGAL::internal_np::maximum_normal_deviation));
  check_same_type<9035>(get_parameter(np, CGAL::internal_np::scan_angle_map));
  check_same_type<9036>(get_parameter(np, CGAL::internal_np::scanline_id_map));
}

int main()
{
  test(CGAL::parameters::vertex_index_map(A<0>(0))
                         .visitor(A<1>(1))
                         .vertex_point_map(A<2>(2))
                         .halfedge_index_map(A<3>(3))
                         .edge_index_map(A<4>(4))
                         .face_index_map(A<5>(5))
                         .edge_is_constrained_map(A<6>(6))
                         .first_index(A<7>(7))
                         .number_of_iterations(A<8>(8))
                         .METIS_options(A<800000001>(800000001))
                         .vertex_partition_id_map(A<800000002>(800000002))
                         .face_partition_id_map(A<800000003>(800000003))
                         .vertex_to_vertex_output_iterator(A<800000004>(800000004))
                         .halfedge_to_halfedge_output_iterator(A<800000005>(800000005))
                         .face_to_face_output_iterator(A<800000006>(800000006))
                         .vertex_to_vertex_map(A<800000007>(800000007))
                         .halfedge_to_halfedge_map(A<800000008>(800000008))
                         .face_to_face_map(A<800000009>(800000009))
                         .implementation_tag(A<800000010>(800000010))
                         .prevent_unselection(A<800000011>(800000011))
                         .stream_precision(A<800000012>(800000012))
                         .vertex_feature_degree_map(A<9>(9))
                         .geom_traits(A<10>(10))
                         .vertex_incident_patches_map(A<11>(11))
                         .density_control_factor(A<12>(12))
                         .use_delaunay_triangulation(A<13>(13))
                         .use_2d_constrained_delaunay_triangulation(A<4573>(4573))
                         .fairing_continuity(A<14>(14))
                         .sparse_linear_solver(A<15>(15))
                         .number_of_relaxation_steps(A<16>(16))
                         .protect_constraints(A<17>(17))
                         .relax_constraints(A<18>(18))
                         .collapse_constraints(A<43>(43))
                         .vertex_is_constrained_map(A<19>(19))
                         .face_patch_map(A<20>(20))
                         .use_random_uniform_sampling(A<21>(21))
                         .use_grid_sampling(A<22>(22))
                         .use_monte_carlo_sampling(A<23>(23))
                         .do_sample_edges(A<24>(24))
                         .do_sample_vertices(A<25>(25))
                         .do_sample_faces(A<26>(26))
                         .number_of_points_on_faces(A<27>(27))
                         .number_of_points_per_face(A<28>(28))
                         .grid_spacing(A<29>(29))
                         .number_of_points_per_edge(A<30>(30))
                         .number_of_points_on_edges(A<31>(31))
                         .number_of_points_per_area_unit(A<32>(32))
                         .number_of_points_per_distance_unit(A<33>(33))
                         .get_cost(A<34>(34))
                         .get_placement(A<35>(35))
                         .face_normal_map(A<36>(36))
                         .random_seed(A<37>(37))
                         .do_project(A<38>(38))
                         .algorithm(A<456>(456))
                         .weight_calculator(A<39>(39))
                         .preserve_genus(A<40>(40))
                         .verbosity_level(A<41>(41))
                         .projection_functor(A<42>(42))
                         .throw_on_self_intersection(A<43>(43))
                         .clip_volume(A<44>(44))
                         .use_compact_clipper(A<45>(45))
                         .non_manifold_feature_map(A<60>(60))
                         .filter(A<61>(61))
                         .face_epsilon_map(A<62>(62))
                         .apply_per_connected_component(A<46>(46))
                         .output_iterator(A<47>(47))
                         .erase_all_duplicates(A<48>(48))
                         .require_same_orientation(A<49>(49))
                         .use_bool_op_to_clip_surface(A<50>(50))
                         .use_binary_mode(A<51>(51))
                         .face_size_map(A<52>(52))
                         .use_angle_smoothing(A<53>(53))
                         .use_area_smoothing(A<54>(54))
                         .use_Delaunay_flips(A<55>(55))
                         .use_safety_constraints(A<56>(56))
                         .do_self_intersection_tests(A<12340>(12340))
                         .do_orientation_tests(A<12341>(12341))
                         .error_codes(A<12342>(12342))
                         .volume_inclusions(A<12343>(12343))
                         .face_connected_component_map(A<12344>(12344))
                         .connected_component_id_to_volume_id(A<12345>(12345))
                         .is_cc_outward_oriented(A<12346>(12346))
                         .intersecting_volume_pairs_output_iterator(A<12347>(12347))
                         .i_used_as_a_predicate(A<12348>(12348))
                         .nesting_levels(A<12349>(12349))
                         .i_used_for_volume_orientation(A<12350>(12350))
                         .area_threshold(A<57>(57))
                         .volume_threshold(A<58>(58))
                         .snapping_tolerance(A<59>(59))
                         .dry_run(A<60>(60))
                         .do_lock_mesh(A<61>(61))
                         .halfedges_keeper(A<62>(62))
                         .use_convex_hull(A<63>(63))
                         .do_simplify_border(A<64>(64))
                         .do_not_modify(A<65>(65))
                         .allow_self_intersections(A<66>(66))
                         .polyhedral_envelope_epsilon(A<67>(67))
                         .point_map(A<9000>(9000))
                         .query_point_map(A<9001>(9001))
                         .normal_map(A<9002>(9002))
                         .diagonalize_traits(A<9003>(9003))
                         .svd_traits(A<9004>(9004))
                         .callback(A<9005>(9005))
                         .sharpness_angle(A<9006>(9006))
                         .edge_sensitivity(A<9007>(9007))
                         .neighbor_radius(A<9008>(9008))
                         .number_of_output_points(A<9009>(9009))
                         .size(A<9010>(9010))
                         .maximum_variation(A<9011>(9011))
                         .degree_fitting(A<9012>(9012))
                         .degree_monge(A<9013>(9013))
                         .threshold_percent(A<9014>(9014))
                         .threshold_distance(A<9015>(9015))
                         .attraction_factor(A<9016>(9016))
                         .plane_map(A<9017>(9017))
                         .plane_index_map(A<9018>(9018))
                         .select_percentage(A<9019>(9019))
                         .require_uniform_sampling(A<9020>(9020))
                         .point_is_constrained_map(A<9021>(9021))
                         .number_of_samples(A<9022>(9022))
                         .accuracy(A<9023>(9023))
                         .maximum_running_time(A<9024>(9024))
                         .overlap(A<9025>(9025))
                         .transformation(A<9026>(9026))
                         .point_set_filters(A<9027>(9027))
                         .matcher(A<9028>(9028))
                         .outlier_filters(A<9029>(9029))
                         .error_minimizer(A<9030>(9030))
                         .transformation_checkers(A<9031>(9031))
                         .inspector(A<9032>(9032))
                         .logger(A<9033>(9033))
                         .maximum_normal_deviation(A<9034>(9034))
                         .scan_angle_map(A<9035>(9035))
                         .scanline_id_map(A<9036>(9036))
                         .maximum_number_of_faces(A<78910>(78910))
       );
  return EXIT_SUCCESS;
}
