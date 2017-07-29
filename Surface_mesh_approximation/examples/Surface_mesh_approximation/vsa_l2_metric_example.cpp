#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;

typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;
typedef CGAL::L2Metric<Polyhedron_3, FacetAreaMap> L2Metric;
typedef CGAL::L2ProxyFitting<Polyhedron_3> L2ProxyFitting;
typedef CGAL::PCAPlaneFitting<Polyhedron_3> PCAPlaneFitting;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return 0;
  
  // read Polyhedron
  Polyhedron_3 mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // facet area map
  std::map<Facet_handle, FT> facet_areas;
  for (Facet_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    Polyhedron_3::Halfedge_handle he = fitr->halfedge();
    const Point_3 &p0 = he->opposite()->vertex()->point();
    const Point_3 &p1 = he->vertex()->point();
    const Point_3 &p2 = he->next()->vertex()->point();

    FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, farea));
  }
  FacetAreaMap area_pmap(facet_areas);

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Point_3> anchor_pos;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;

  CGAL::vsa_extract(mesh, tris, anchor_pos,
    PCAPlaneFitting(mesh),
    L2Metric(mesh, area_pmap),
    L2ProxyFitting(mesh),
    init, num_proxies, num_iterations);

  return EXIT_SUCCESS;
}
