#include <cassert>
#include "test_meshing_utilities.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>

#include <sstream>

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Polyhedron_tester : public Tester<K>
{
  void polyhedron() const
  {
    // Domain
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
    typedef typename CGAL::Mesh_polyhedron_3<K>::type MeshPolyhedron_3;

    // Triangulation
    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,
      typename Mesh_domain::Corner_index,
      typename Mesh_domain::Curve_index> C3t3;

    // Criteria
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    MeshPolyhedron_3 P;

    std::stringstream poly("OFF\n"
                           "4 4 0\n"
                           "0 0 0\n"
                           "0 0 1\n"
                           "0 1 0\n"
                           "1 0 0\n"
                           "3  3 1 2\n"
                           "3  0 1 3\n"
                           "3  0 3 2\n"
                           "3  0 2 1\n");
    poly >> P;
    assert(P.size_of_vertices() == 4);
    assert(P.size_of_facets() == 4);
    Mesh_domain domain(P);

    // Get sharp features
    domain.detect_features();

    // Mesh criteria
    using namespace CGAL::parameters;
    const double cs = 0.408248;
    Mesh_criteria criteria(
#ifndef CGAL_MESH_3_USE_DEFAULT_EDGE_SIZE
                           edge_size = cs/2.0,
#endif // see test_meshing_with_default_edge_size.cpp
                           facet_angle = 25,
                           facet_size = cs,
                           facet_distance = cs/10.0,
                           cell_radius_edge_ratio = 3,
                           cell_size = cs);
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    CGAL::remove_far_points_in_mesh_3(c3t3);

    double vol = 1/6.;
    this->verify_c3t3_volume(c3t3, vol*0.95, vol*1.05);
#ifdef CGAL_MESH_3_USE_DEFAULT_EDGE_SIZE
    this->verify_c3t3(c3t3,domain,Polyhedral_tag(),
                      55, 70, 110, 135, 85, 130);
#else // not CGAL_MESH_3_USE_DEFAULT_EDGE_SIZE
    this->verify_c3t3(c3t3,domain,Polyhedral_tag(),
                      55, 65, 110, 125, 85, 120);
#endif // not CGAL_MESH_3_USE_DEFAULT_EDGE_SIZE

#ifndef DUMP_FILES_PREFIX
#  define DUMP_FILES_PREFIX "unit_tetrahedron"
#endif// see test_meshing_with_default_edge_size.cpp
    CGAL::dump_c3t3(c3t3, DUMP_FILES_PREFIX);
  }
};

int main() {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  Polyhedron_tester<K> test_epic;
  std::cerr << "Mesh generation from a polyhedron:\n";
  test_epic.polyhedron();

#ifdef CGAL_LINKED_WITH_TBB
  Polyhedron_tester<K, CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a polyhedron:\n";
  test_epic_p.polyhedron();
#endif

  return EXIT_SUCCESS;
}
