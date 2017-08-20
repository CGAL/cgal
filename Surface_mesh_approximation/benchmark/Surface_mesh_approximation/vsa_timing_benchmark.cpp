#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/Timer.h>

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

typedef CGAL::Timer Timer;

/**
 * This file is a timing benchmark of each phase.
 * With different configuration:
   1. initialization
   2. number of proxies
   3. number of iterations
 * TODO: error
 */
int main(int argc, char *argv[])
{
  if (argc < 5)
    return 1;

  Polyhedron_3 mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return 1;
  }
  std::cerr << "#triangles " << mesh.size_of_facets() << std::endl;

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

  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return 1;
  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::cerr << "#init " << init << std::endl;
  std::cerr << "#num_proxies " << num_proxies << std::endl;
  std::cerr << "#num_iterations " << num_iterations << std::endl;

  Timer t0, t1;
  t1.start();

  std::cerr << "start initialization" << std::endl;
  t0.reset();
  t0.start();
  vsa_l21.init_proxies(num_proxies, static_cast<VSAL21::Initialization>(init));
  t0.stop();
  std::cerr << "initialization time " << t0.time() << " sec." << std::endl;

  std::cerr << "start relaxation" << std::endl;
  t0.reset();
  t0.start();
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();
  t0.stop();
  std::cerr << "relaxation time " << t0.time() << " sec." << std::endl;

  Polyhedron_3 mesh_out;
  std::cerr << "start relaxation" << std::endl;
  t0.reset();
  t0.start();
  vsa_l21.meshing(mesh_out);
  t0.stop();
  std::cerr << "meshing time " << t0.time() << " sec." << std::endl;

  std::cerr << "total time " << t1.time() << " sec." << std::endl;

  return 0;
}
