#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Detect_features_in_image.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Sizing_field_with_aabb_tree.h>

#include "test_meshing_utilities.h"

#define CGAL_MESH_3_VERBOSE true

struct Edge_distance_test_helper
{
  typedef K_e_i K;
  typedef K_e_i::FT FT;

  template <typename C3t3, typename Mesh_Domain>
  std::tuple<FT, FT, FT> compute_stats(const C3t3 & c3t3, const Mesh_Domain & domain) {
    typename C3t3::Edges_in_complex_iterator edge_begin = c3t3.edges_in_complex_begin();
    FT min_edge_size = CGAL::sqrt(
      CGAL::squared_distance(edge_begin->first->vertex( edge_begin->second )->point(), edge_begin->first->vertex( edge_begin->third )->point())
      );
    FT avg_edge_size = 0;
    FT sum_approx_error = 0;
    int nbEdges = 0;
    for (typename C3t3::Edges_in_complex_iterator eit = edge_begin ; eit != c3t3.edges_in_complex_end(); ++eit ) {
      const typename C3t3::Vertex_handle& va = eit->first->vertex(eit->second);
      const typename C3t3::Vertex_handle& vb = eit->first->vertex(eit->third);

      // Get edge distance

      FT dist = CGAL::sqrt(
        CGAL::squared_distance(va->point(), vb->point())
      );
      if (dist < min_edge_size) {
        min_edge_size = dist;
      }
      avg_edge_size += dist;
      nbEdges++;

      // Get edge approximation error

      typename C3t3::Curve_index curve_index = domain.curve_index((va->in_dimension() < vb->in_dimension()) ? vb->index() : va->index());

      const K::Point_3& pa = va->point().point();
      const K::Point_3& pb = vb->point().point();
      const K::Point_3& segment_middle = CGAL::midpoint(pa, pb);
      // Obtain the geodesic middle point
      FT signed_geodesic_distance = domain.signed_geodesic_distance(pa, pb, curve_index);
      K::Point_3 geodesic_middle;
      if (signed_geodesic_distance >= 0)
      {
        geodesic_middle = domain.construct_point_on_curve(pa, curve_index, signed_geodesic_distance / 2);
      }
      else
      {
        geodesic_middle = domain.construct_point_on_curve(pb, curve_index, -signed_geodesic_distance / 2);
      }
      // Compare distance to the parameter's distance
      sum_approx_error += CGAL::squared_distance(segment_middle, geodesic_middle);
    }
    avg_edge_size /= nbEdges;

    return std::make_tuple(min_edge_size, avg_edge_size, sum_approx_error);
  }

  /**
  * @brief verify that there are more vertices than without criteria
  * and that the edges are smaller (minimum size and average size)
  * and that the approximation of polylines features is better
  */
  template <typename C3t3, typename Mesh_Domain>
  void test_c3t3_without_and_with(const C3t3 & c3t3_without, const Mesh_Domain & domain_without, const C3t3 & c3t3_with, const Mesh_Domain & domain_with) {
    std::tuple<FT, FT, FT> stats_without = this->compute_stats(c3t3_without, domain_without);
    FT min_edge_size_without_distance = std::get<0>(stats_without);
    FT avg_edge_size_without_distance = std::get<1>(stats_without);
    FT sum_approx_error_without_distance = std::get<2>(stats_without);
    std::size_t nb_vertices_without_distance  = c3t3_without.triangulation().number_of_vertices();
    std::size_t nb_triangles_without_distance = c3t3_without.number_of_facets_in_complex();

    std::tuple<FT, FT, FT> stats_with= this->compute_stats(c3t3_with, domain_with);
    FT min_edge_size_with_distance = std::get<0>(stats_with);
    FT avg_edge_size_with_distance = std::get<1>(stats_with);
    FT sum_approx_error_with_distance = std::get<2>(stats_with);
    std::size_t nb_vertices_with_distance  = c3t3_with.triangulation().number_of_vertices();
    std::size_t nb_triangles_with_distance = c3t3_with.number_of_facets_in_complex();

    assert(nb_vertices_with_distance  > nb_vertices_without_distance);
    assert(nb_triangles_with_distance > nb_triangles_without_distance);

    assert(min_edge_size_with_distance < min_edge_size_without_distance);
    assert(avg_edge_size_with_distance < avg_edge_size_without_distance);

    assert(sum_approx_error_with_distance < sum_approx_error_without_distance);
  }
};

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Distance_polyhedral_tester : public Tester<K_e_i>, public Edge_distance_test_helper
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  // Domain
  typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;

  // Triangulation
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, K, Concurrency_tag>::type Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

public:
  void operator()(const std::size_t min_vertices_expected = 0,
                  const std::size_t max_vertices_expected = STD_SIZE_T_MAX,
                  const std::size_t min_facets_expected = 0,
                  const std::size_t max_facets_expected = STD_SIZE_T_MAX)
  {
    const std::string fname = CGAL::data_file_path("meshes/dragknob.off");

    std::ifstream input(fname);
    using namespace CGAL::parameters;

    Polyhedron polyhedron;
    input >> polyhedron;
    if (input.fail()) {
      std::cerr << "Error: Cannot read file " << fname << std::endl;
      return;
    }

    if (!CGAL::is_triangle_mesh(polyhedron)) {
      std::cerr << "Input geometry is not triangulated." << std::endl;
      return;
    }

    // Create domain
    Mesh_domain domain(polyhedron);

    domain.detect_features(40);

    // Mesh criteria
    Mesh_criteria criteria(edge_size = 0.074,
        edge_distance = 0.00074,
        facet_distance = 0.0074,
        facet_angle = 25,
        facet_size = 0.074,
        cell_radius_edge_ratio = 3,
        cell_size = 0.074);

    assert(criteria.edge_criteria_object().has_distance_field());

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    this->verify_c3t3(c3t3, domain, Polyhedral_tag(),
                      min_vertices_expected, max_vertices_expected,
                      min_facets_expected, max_facets_expected);

    Mesh_criteria criteria_without(edge_size = 0.074,
        facet_distance = 0.0074,
        facet_angle = 25,
        facet_size = 0.074,
        cell_radius_edge_ratio = 3,
        cell_size = 0.074);

    assert(!criteria_without.edge_criteria_object().has_distance_field());

    // Mesh generation
    C3t3 c3t3_without = CGAL::make_mesh_3<C3t3>(domain, criteria_without, no_perturb(), no_exude());

    this->test_c3t3_without_and_with(c3t3_without, domain, c3t3, domain);
  }
};

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Distance_label_image_tester : public Tester<K_e_i>, public Edge_distance_test_helper
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  // Domain
  typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
  typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;

  // Triangulation
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, K, Concurrency_tag>::type Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

public:

  void operator()(const std::size_t min_vertices_expected = 0,
                  const std::size_t max_vertices_expected = STD_SIZE_T_MAX,
                  const std::size_t min_facets_expected = 0,
                  const std::size_t max_facets_expected = STD_SIZE_T_MAX)
  {
    const std::string fname = CGAL::data_file_path("images/quadDomainCube.inr");

    using namespace CGAL::parameters;

    CGAL::Image_3 image;
    if(!image.read(fname)){
        std::cerr << "Error: Cannot read file " <<  fname << std::endl;
        return;
    }

    // Create domain
    Mesh_domain domain
        = Mesh_domain::create_labeled_image_mesh_domain(image,
                features_detector = CGAL::Mesh_3::Detect_features_in_image());

    // Mesh criteria
    Mesh_criteria criteria(edge_size = 5.,
        edge_distance = 0.3,
        facet_distance = 0.3,
        facet_angle = 25,
        facet_size = 5.,
        cell_radius_edge_ratio = 3,
        cell_size = 5.);

    assert(criteria.edge_criteria_object().has_distance_field());

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    this->verify_c3t3(c3t3, domain, Polyhedral_tag(),
                      min_vertices_expected, max_vertices_expected,
                      min_facets_expected, max_facets_expected);

    // verify that there are more vertices than without criteria
    Mesh_criteria criteria_without(edge_size = 5.,
        facet_distance = 0.3,
        facet_angle = 25,
        facet_size = 5.,
        cell_radius_edge_ratio = 3,
        cell_size = 5.);

    assert(!criteria_without.edge_criteria_object().has_distance_field());

    // Mesh generation
    C3t3 c3t3_without = CGAL::make_mesh_3<C3t3>(domain, criteria_without, no_perturb(), no_exude());

    this->test_c3t3_without_and_with(c3t3_without, domain, c3t3, domain);
  }
};


int main()
{
  Distance_polyhedral_tester<> test_epic_poly;
  std::cerr << "Mesh generation with edge_distance from polyhedral domain:\n";
  test_epic_poly(1820 * 0.95, 1820 * 1.05, 3012 * 0.95, 3012 * 1.05);

  Distance_label_image_tester<> test_epic_image;
  std::cerr << "Mesh generation with edge_distance from label image domain:\n";
  test_epic_image(776 * 0.95, 776 * 1.05, 1525 * 0.95, 1525 * 1.05);

#ifdef CGAL_LINKED_WITH_TBB
  Distance_polyhedral_tester<CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation with edge_distance from polyhedral domain:\n";
  test_epic_p();

  Distance_label_image_tester<CGAL::Parallel_tag> test_epic_image_p;
  std::cerr << "Parallel mesh generation with edge_distance from label image domain:\n";
  test_epic_image_p();
#endif
  return EXIT_SUCCESS;
}
