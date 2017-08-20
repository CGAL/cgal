#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_handle, Vector_3> > FacetNormalMap;

typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;
typedef CGAL::L21Metric<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21Metric;
typedef CGAL::L21ProxyFitting<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L21Metric, L21ProxyFitting> VSAL21;

bool test_shape(const std::string file_name, const std::size_t target_num_proxies)
{
  Polyhedron_3 mesh;
  std::ifstream input(file_name);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return false;
  }

  std::map<Facet_handle, Vector_3> facet_normals;
  std::map<Facet_handle, FT> facet_areas;
  for (Polyhedron_3::Facet_iterator fitr = mesh.facets_begin();
    fitr != mesh.facets_end(); ++fitr) {
    Polyhedron_3::Halfedge_handle he = fitr->halfedge();
    const Point_3 &p0 = he->opposite()->vertex()->point();
    const Point_3 &p1 = he->vertex()->point();
    const Point_3 &p2 = he->next()->vertex()->point();

    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<Facet_handle, Vector_3>(fitr, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  L21Metric l21_metric(normal_pmap, area_pmap);
  L21ProxyFitting l21_fitting(normal_pmap, area_pmap);

  // algorithm instance
  VSAL21 vsa_l21(l21_metric, l21_fitting);
  vsa_l21.set_mesh(mesh);

  // approximation, init from error, drop to the target error incrementally
  // should reach targeted number of proxies gradually
  const FT drop(1e-8);
  const std::size_t num_iterations = 20;
  vsa_l21.init_proxies_error(drop, VSAL21::IncrementalInit);
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();
  if (vsa_l21.get_proxies_size() != target_num_proxies) {
    std::cout << "#targeted - #result "
      << target_num_proxies << ' '
      << vsa_l21.get_proxies_size() << std::endl;
    std::cout << "incremental reaching failed" << std::endl;
    return false;
  }

  // meshing, should be manifold
  Polyhedron_3 mesh_out;
  if (!vsa_l21.meshing(mesh_out)) {
    std::cout << "incremental reaching meshing non-manifold" << std::endl;
    return false;
  }

  return true;
}

/**
 * This file tests the correctness of the algorithm.
 * Basically we input a cube mesh and see if it outputs a cube.
 */
int main()
{
  const std::string file0("./data/cube_meshed.off");
  std::cout << "Testing file \"" << file0 << '\"';
  if (!test_shape(file0, 6)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  const std::string file1("./data/cube_meshed_open.off");
  std::cout << "Testing file \"" << file1 << '\"';
  if (!test_shape(file1, 5)) {
    std::cout << "Failed." << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "Succeeded." << std::endl;

  return EXIT_SUCCESS;
}
