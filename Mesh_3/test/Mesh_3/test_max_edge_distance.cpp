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

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Distance_polyhedral_tester : public Tester<K_e_i>
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
  void operator()(std::size_t expected_nb_vertices, std::size_t expected_nb_triangles)
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

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    std::size_t nb_vertices  = c3t3.triangulation().number_of_vertices();
    std::size_t nb_triangles = c3t3.number_of_facets_in_complex();

    this->verify_c3t3(c3t3, domain, Polyhedral_tag(), expected_nb_vertices * 0.95, expected_nb_vertices * 1.05, expected_nb_triangles * 0.95, expected_nb_triangles * 1.05);

    // verify that there are more vertices than without criteria
    // and that the edges are smaller (minimum size and average size)

    C3t3::Edges_in_complex_iterator edge_begin = c3t3.edges_in_complex_begin();
    K::FT min_edge_size_with_distance = CGAL::sqrt(
        CGAL::squared_distance(edge_begin->first->vertex( edge_begin->second )->point(), edge_begin->first->vertex( edge_begin->third )->point())
        );
    K::FT avg_edge_size_with_distance = min_edge_size_with_distance;
    edge_begin++;
    int nbEdges = 1;
    for (C3t3::Edges_in_complex_iterator eit = edge_begin ; eit != c3t3.edges_in_complex_end(); ++eit ) {
      K::FT dist = CGAL::sqrt(
          CGAL::squared_distance(eit->first->vertex( eit->second )->point(), eit->first->vertex( eit->third )->point())
          );
      if (dist < min_edge_size_with_distance) {
          min_edge_size_with_distance = dist;
      }
      avg_edge_size_with_distance += dist;
      nbEdges++;
    }
    avg_edge_size_with_distance /= nbEdges;

    Mesh_criteria criteria_without(edge_size = 0.074,
        facet_distance = 0.0074,
        facet_angle = 25,
        facet_size = 0.074,
        cell_radius_edge_ratio = 3,
        cell_size = 0.074);

    // Mesh generation
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria_without, no_perturb(), no_exude());

    // Get edge min and average length
    C3t3::Edges_in_complex_iterator edge_begin_ = c3t3.edges_in_complex_begin();
    K::FT min_edge_size_without_distance = CGAL::sqrt(
        CGAL::squared_distance(edge_begin_->first->vertex( edge_begin_->second )->point(), edge_begin_->first->vertex( edge_begin_->third )->point())
        );
    K::FT avg_edge_size_without_distance = min_edge_size_without_distance;
    edge_begin_++;
    nbEdges = 1;
    for (C3t3::Edges_in_complex_iterator eit = edge_begin_ ; eit != c3t3.edges_in_complex_end(); ++eit ) {
      K::FT dist = CGAL::sqrt(
          CGAL::squared_distance(eit->first->vertex( eit->second )->point(), eit->first->vertex( eit->third )->point())
          );
      if (dist < min_edge_size_without_distance) {
          min_edge_size_without_distance = dist;
      }
      avg_edge_size_without_distance += dist;
      nbEdges++;
    }
    avg_edge_size_without_distance /= nbEdges;

    assert(nb_vertices  > c3t3.triangulation().number_of_vertices());
    assert(nb_triangles > c3t3.number_of_facets_in_complex());

    assert(min_edge_size_with_distance < min_edge_size_without_distance);
    assert(avg_edge_size_with_distance < avg_edge_size_without_distance);
  }
};

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Distance_label_image_tester : public Tester<K_e_i>
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

  void operator()(std::size_t expected_nb_vertices, std::size_t expected_nb_triangles)
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

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    std::size_t nb_vertices  = c3t3.triangulation().number_of_vertices();
    std::size_t nb_triangles = c3t3.number_of_facets_in_complex();

    this->verify_c3t3(c3t3, domain, Polyhedral_tag(), expected_nb_vertices * 0.95, expected_nb_vertices * 1.05, expected_nb_triangles * 0.95, expected_nb_triangles * 1.05);

    // verify that there are more vertices than without criteria
    Mesh_criteria criteria_without(edge_size = 5.,
        facet_distance = 0.3,
        facet_angle = 25,
        facet_size = 5.,
        cell_radius_edge_ratio = 3,
        cell_size = 5.);

    // Mesh generation
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria_without, no_perturb(), no_exude());
    assert(nb_vertices  > c3t3.triangulation().number_of_vertices());
    assert(nb_triangles > c3t3.number_of_facets_in_complex());
  }
};


int main(int argc, char* argv[])
{
  Distance_polyhedral_tester<> test_epic_poly;
  std::cerr << "Mesh generation with edge_distance from polyhedral domain:\n";
  test_epic_poly(1820, 3012);

  Distance_label_image_tester<> test_epic_image;
  std::cerr << "Mesh generation with edge_distance from label image domain:\n";
  test_epic_image(776, 1525);

#ifdef CGAL_LINKED_WITH_TBB
  Distance_polyhedral_tester<CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation with edge_distance from polyhedral domain:\n";
  test_epic_p(1924, 3034);

  Distance_label_image_tester<CGAL::Parallel_tag> test_epic_image_p;
  std::cerr << "Parallel mesh generation with edge_distance from label image domain:\n";
  test_epic_image_p(847, 1509);
#endif
  return EXIT_SUCCESS;
}
