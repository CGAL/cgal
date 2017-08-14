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
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;

typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;
typedef boost::associative_property_map<std::map<Facet_handle, Vector_3> > FacetNormalMap;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef CGAL::L21Metric<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21Metric;
typedef CGAL::L21ProxyFitting<Polyhedron_3, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L21Metric, L21ProxyFitting> VSA;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return EXIT_FAILURE;

  // create and read Polyhedron_3
  Polyhedron_3 mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;

  // construct facet normal & area map
  std::map<Facet_handle, Vector_3> facet_normals;
  std::map<Facet_handle, FT> facet_areas;
  for (Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Polyhedron_3::Halfedge_around_facet_circulator he = fitr->facet_begin();
    const Point_3 p0 = he->opposite()->vertex()->point();
    const Point_3 p1 = he->vertex()->point();
    const Point_3 p2 = he->next()->vertex()->point();
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<Facet_handle, Vector_3>(fitr, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  VSA vsa_approx(L21Metric(normal_pmap, area_pmap), L21ProxyFitting(normal_pmap, area_pmap));

  vsa_approx.set_mesh(mesh);
  vsa_approx.init_proxies(num_proxies, VSA::RandomInit);
  
  for (std::size_t i = 0; i < num_iterations; ++i) {
    vsa_approx.partition();
    vsa_approx.fit();
  }

  return EXIT_SUCCESS;
}
