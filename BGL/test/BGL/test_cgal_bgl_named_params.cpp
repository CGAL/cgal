#include <CGAL/boost/graph/named_function_params.h>
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
  using boost::get_param;

  // Test values

    // Named parameters from Boost
  assert(get_param(np, boost::vertex_index).v == 0);
  assert(get_param(np, boost::graph_visitor).v == 1);

    // Named parameters that we use in CGAL
  assert(get_param(np, CGAL::internal_np::vertex_point).v == 2);
  assert(get_param(np, CGAL::internal_np::halfedge_index).v == 3);
  assert(get_param(np, CGAL::internal_np::edge_index).v == 4);
  assert(get_param(np, CGAL::internal_np::face_index).v == 5);

  assert(get_param(np, CGAL::internal_np::edge_is_constrained).v == 6);
  assert(get_param(np, CGAL::internal_np::first_index).v == 7);
  assert(get_param(np, CGAL::internal_np::number_of_iterations).v == 8);

  assert(get_param(np, CGAL::internal_np::METIS_options).v == 800000001);
  assert(get_param(np, CGAL::internal_np::vertex_partition_id).v == 800000002);
  assert(get_param(np, CGAL::internal_np::face_partition_id).v == 800000003);
  
  assert(get_param(np, CGAL::internal_np::vertex_to_vertex_output_iterator).v == 800000004);
  assert(get_param(np, CGAL::internal_np::halfedge_to_halfedge_output_iterator).v == 800000005);
  assert(get_param(np, CGAL::internal_np::face_to_face_output_iterator).v == 800000006);
  
  assert(get_param(np, CGAL::internal_np::vertex_to_vertex_map).v == 800000007);
  assert(get_param(np, CGAL::internal_np::halfedge_to_halfedge_map).v == 800000008);
  assert(get_param(np, CGAL::internal_np::face_to_face_map).v == 800000009);
 

    // Named parameters that we use in the package 'Mesh_3'
  assert(get_param(np, CGAL::internal_np::vertex_feature_degree).v == 9);

    // Named parameters used in the package 'Polygon Mesh Processing'
  assert(get_param(np, CGAL::internal_np::geom_traits).v == 10);
  assert(get_param(np, CGAL::internal_np::vertex_incident_patches).v == 11);
  assert(get_param(np, CGAL::internal_np::density_control_factor).v == 12);
  assert(get_param(np, CGAL::internal_np::use_delaunay_triangulation).v == 13);
  assert(get_param(np, CGAL::internal_np::fairing_continuity).v == 14);
  assert(get_param(np, CGAL::internal_np::sparse_linear_solver).v == 15);
  assert(get_param(np, CGAL::internal_np::number_of_relaxation_steps).v == 16);
  assert(get_param(np, CGAL::internal_np::protect_constraints).v == 17);
  assert(get_param(np, CGAL::internal_np::relax_constraints).v == 18);
  assert(get_param(np, CGAL::internal_np::collapse_constraints).v == 43);
  assert(get_param(np, CGAL::internal_np::vertex_is_constrained).v == 19);
  assert(get_param(np, CGAL::internal_np::face_patch).v == 20);
  assert(get_param(np, CGAL::internal_np::random_uniform_sampling).v == 21);
  assert(get_param(np, CGAL::internal_np::grid_sampling).v == 22);
  assert(get_param(np, CGAL::internal_np::monte_carlo_sampling).v == 23);
  assert(get_param(np, CGAL::internal_np::do_sample_edges).v == 24);
  assert(get_param(np, CGAL::internal_np::do_sample_vertices).v == 25);
  assert(get_param(np, CGAL::internal_np::do_sample_faces).v == 26);
  assert(get_param(np, CGAL::internal_np::number_of_points_on_faces).v == 27);
  assert(get_param(np, CGAL::internal_np::number_of_points_per_face).v == 28);
  assert(get_param(np, CGAL::internal_np::grid_spacing).v == 29);
  assert(get_param(np, CGAL::internal_np::number_of_points_per_edge).v == 30);
  assert(get_param(np, CGAL::internal_np::number_of_points_on_edges).v == 31);
  assert(get_param(np, CGAL::internal_np::nb_points_per_area_unit).v == 32);
  assert(get_param(np, CGAL::internal_np::nb_points_per_distance_unit).v == 33);
  assert(get_param(np, CGAL::internal_np::throw_on_self_intersection).v == 43);
  assert(get_param(np, CGAL::internal_np::clip_volume).v == 44);
  assert(get_param(np, CGAL::internal_np::use_compact_clipper).v == 45);
  assert(get_param(np, CGAL::internal_np::erase_all_duplicates).v == 48);
  assert(get_param(np, CGAL::internal_np::require_same_orientation).v == 49);
  assert(get_param(np, CGAL::internal_np::use_bool_op_to_clip_surface).v == 50);
  assert(get_param(np, CGAL::internal_np::gradient_descent_precision).v == 54);
  assert(get_param(np, CGAL::internal_np::use_explicit_scheme).v == 55);


    // Named parameters that we use in the package 'Surface Mesh Simplification'
  assert(get_param(np, CGAL::internal_np::get_cost_policy).v == 34);
  assert(get_param(np, CGAL::internal_np::get_placement_policy).v == 35);

    // To-be-documented named parameters
  assert(get_param(np, CGAL::internal_np::face_normal).v == 36);
  assert(get_param(np, CGAL::internal_np::random_seed).v == 37);
  assert(get_param(np, CGAL::internal_np::do_project).v == 38);

    // Internal named parameters
  assert(get_param(np, CGAL::internal_np::weight_calculator).v == 39);
  assert(get_param(np, CGAL::internal_np::preserve_genus).v == 40);
  assert(get_param(np, CGAL::internal_np::verbosity_level).v == 41);
  assert(get_param(np, CGAL::internal_np::use_binary_mode).v == 51);
  assert(get_param(np, CGAL::internal_np::projection_functor).v == 42);
  assert(get_param(np, CGAL::internal_np::apply_per_connected_component).v == 46);
  assert(get_param(np, CGAL::internal_np::output_iterator).v == 47);


  // Test types

    // Named parameters from Boost
  check_same_type<0>(get_param(np, boost::vertex_index));
  check_same_type<1>(get_param(np, boost::graph_visitor));

    // Named parameters that we use in CGAL
  check_same_type<2>(get_param(np, CGAL::internal_np::vertex_point));
  check_same_type<3>(get_param(np, CGAL::internal_np::halfedge_index));
  check_same_type<4>(get_param(np, CGAL::internal_np::edge_index));
  check_same_type<5>(get_param(np, CGAL::internal_np::face_index));

  check_same_type<6>(get_param(np, CGAL::internal_np::edge_is_constrained));
  check_same_type<7>(get_param(np, CGAL::internal_np::first_index));
  check_same_type<8>(get_param(np, CGAL::internal_np::number_of_iterations));

  check_same_type<800000001>(get_param(np, CGAL::internal_np::METIS_options));
  check_same_type<800000002>(get_param(np, CGAL::internal_np::vertex_partition_id));
  check_same_type<800000003>(get_param(np, CGAL::internal_np::face_partition_id));
  check_same_type<800000004>(get_param(np, CGAL::internal_np::vertex_to_vertex_output_iterator));
  check_same_type<800000005>(get_param(np, CGAL::internal_np::halfedge_to_halfedge_output_iterator));
  check_same_type<800000006>(get_param(np, CGAL::internal_np::face_to_face_output_iterator));
  check_same_type<800000007>(get_param(np, CGAL::internal_np::vertex_to_vertex_map));
  check_same_type<800000008>(get_param(np, CGAL::internal_np::halfedge_to_halfedge_map));
  check_same_type<800000009>(get_param(np, CGAL::internal_np::face_to_face_map));

    // Named parameters that we use in the package 'Mesh_3'
  check_same_type<9>(get_param(np, CGAL::internal_np::vertex_feature_degree));

    // Named parameters used in the package 'Polygon Mesh Processing'
  check_same_type<10>(get_param(np, CGAL::internal_np::geom_traits));
  check_same_type<11>(get_param(np, CGAL::internal_np::vertex_incident_patches));
  check_same_type<12>(get_param(np, CGAL::internal_np::density_control_factor));
  check_same_type<13>(get_param(np, CGAL::internal_np::use_delaunay_triangulation));
  check_same_type<14>(get_param(np, CGAL::internal_np::fairing_continuity));
  check_same_type<15>(get_param(np, CGAL::internal_np::sparse_linear_solver));
  check_same_type<16>(get_param(np, CGAL::internal_np::number_of_relaxation_steps));
  check_same_type<17>(get_param(np, CGAL::internal_np::protect_constraints));
  check_same_type<18>(get_param(np, CGAL::internal_np::relax_constraints));
  check_same_type<43>(get_param(np, CGAL::internal_np::collapse_constraints));
  check_same_type<19>(get_param(np, CGAL::internal_np::vertex_is_constrained));
  check_same_type<20>(get_param(np, CGAL::internal_np::face_patch));
  check_same_type<21>(get_param(np, CGAL::internal_np::random_uniform_sampling));
  check_same_type<22>(get_param(np, CGAL::internal_np::grid_sampling));
  check_same_type<23>(get_param(np, CGAL::internal_np::monte_carlo_sampling));
  check_same_type<24>(get_param(np, CGAL::internal_np::do_sample_edges));
  check_same_type<25>(get_param(np, CGAL::internal_np::do_sample_vertices));
  check_same_type<26>(get_param(np, CGAL::internal_np::do_sample_faces));
  check_same_type<27>(get_param(np, CGAL::internal_np::number_of_points_on_faces));
  check_same_type<28>(get_param(np, CGAL::internal_np::number_of_points_per_face));
  check_same_type<29>(get_param(np, CGAL::internal_np::grid_spacing));
  check_same_type<30>(get_param(np, CGAL::internal_np::number_of_points_per_edge));
  check_same_type<31>(get_param(np, CGAL::internal_np::number_of_points_on_edges));
  check_same_type<32>(get_param(np, CGAL::internal_np::nb_points_per_area_unit));
  check_same_type<33>(get_param(np, CGAL::internal_np::nb_points_per_distance_unit));
  check_same_type<43>(get_param(np, CGAL::internal_np::throw_on_self_intersection));
  check_same_type<44>(get_param(np, CGAL::internal_np::clip_volume));
  check_same_type<45>(get_param(np, CGAL::internal_np::use_compact_clipper));
  check_same_type<48>(get_param(np, CGAL::internal_np::erase_all_duplicates));
  check_same_type<49>(get_param(np, CGAL::internal_np::require_same_orientation));
  check_same_type<50>(get_param(np, CGAL::internal_np::use_bool_op_to_clip_surface));
  check_same_type<54>(get_param(np, CGAL::internal_np::gradient_descent_precision));
  check_same_type<55>(get_param(np, CGAL::internal_np::use_explicit_scheme));

    // Named parameters that we use in the package 'Surface Mesh Simplification'
  check_same_type<34>(get_param(np, CGAL::internal_np::get_cost_policy));
  check_same_type<35>(get_param(np, CGAL::internal_np::get_placement_policy));

    // To-be-documented named parameters
  check_same_type<36>(get_param(np, CGAL::internal_np::face_normal));
  check_same_type<37>(get_param(np, CGAL::internal_np::random_seed));
  check_same_type<38>(get_param(np, CGAL::internal_np::do_project));

    // Internal named parameters
  check_same_type<39>(get_param(np, CGAL::internal_np::weight_calculator));
  check_same_type<40>(get_param(np, CGAL::internal_np::preserve_genus));
  check_same_type<41>(get_param(np, CGAL::internal_np::verbosity_level));
  check_same_type<51>(get_param(np, CGAL::internal_np::use_binary_mode));
  check_same_type<42>(get_param(np, CGAL::internal_np::projection_functor));
  check_same_type<46>(get_param(np, CGAL::internal_np::apply_per_connected_component));
  check_same_type<47>(get_param(np, CGAL::internal_np::output_iterator));
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
                         .vertex_feature_degree_map(A<9>(9))
                         .geom_traits(A<10>(10))
                         .vertex_incident_patches_map(A<11>(11))
                         .density_control_factor(A<12>(12))
                         .use_delaunay_triangulation(A<13>(13))
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
                         .weight_calculator(A<39>(39))
                         .preserve_genus(A<40>(40))
                         .verbosity_level(A<41>(41))
                         .use_binary_mode(A<51>(51))
                         .projection_functor(A<42>(42))
                         .throw_on_self_intersection(A<43>(43))
                         .clip_volume(A<44>(44))
                         .use_compact_clipper(A<45>(45))
                         .use_bool_op_to_clip_surface(A<50>(50))
                         .apply_per_connected_component(A<46>(46))
                         .output_iterator(A<47>(47))
                         .erase_all_duplicates(A<48>(48))
                         .gradient_descent_precision(A<54>(54))
                         .use_explicit_scheme(A<55>(55))
                         .require_same_orientation(A<52>(52))
       );

  return EXIT_SUCCESS;
}
