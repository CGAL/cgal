#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>
#include <CGAL/vsa_mesh_approximation.h>
#include <CGAL/vsa_mesh_approximation_traits.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::FT FT;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_const_handle Facet_const_handle;
typedef Polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator Facet_const_iterator;
typedef boost::associative_property_map<std::map<Facet_const_handle, Vector_3> > FacetNormalMap;
typedef boost::associative_property_map<std::map<Facet_const_handle, FT> > FacetAreaMap;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

typedef CGAL::PlaneProxy<Polyhedron> PlaneProxy;
typedef CGAL::L21Metric<PlaneProxy, FacetNormalMap, FacetAreaMap> L21Metric;
typedef CGAL::L21ProxyFitting<PlaneProxy, L21Metric, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
typedef CGAL::PlaneFitting<Polyhedron, FacetAreaMap, FacetNormalMap, VertexPointMap> PlaneFitting;
typedef CGAL::L21ApproximationTrait<PlaneProxy, Polyhedron, L21Metric, L21ProxyFitting, PlaneFitting, VertexPointMap, FacetNormalMap, FacetAreaMap> L21ApproximationTrait;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return 0;

  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // construct facet normal & area map
  std::map<Facet_const_handle, Vector_3> facet_normals;
  std::map<Facet_const_handle, FT> facet_areas;
  for(Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Halfedge_const_handle he = fitr->halfedge();
    const Point_3 p1 = he->opposite()->vertex()->point();
    const Point_3 p2 = he->vertex()->point();
    const Point_3 p3 = he->next()->vertex()->point();
    Vector_3 normal = CGAL::unit_normal(p1, p2, p3);
    facet_normals.insert(std::pair<Facet_const_handle, Vector_3>(fitr, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p1, p2, p3))));
    facet_areas.insert(std::pair<Facet_const_handle, FT>(fitr, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);
  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<Polyhedron &>(mesh));

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;
  CGAL::vsa_mesh_approximation(mesh,
    tris,
    anchor_pos,
    L21ApproximationTrait(mesh, point_pmap, normal_pmap, area_pmap),
    init,
    num_proxies,
    num_iterations);

  return EXIT_SUCCESS;
}
