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
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Facet_const_handle Facet_const_handle;
typedef Polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator Facet_const_iterator;
typedef boost::associative_property_map<std::map<Facet_const_handle, Vector_3> > FacetNormalMap;
typedef boost::associative_property_map<std::map<Facet_const_handle, FT> > FacetAreaMap;
typedef boost::associative_property_map<std::map<Facet_const_handle, Point_3> > FacetCenterMap;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type VertexPointMap;

struct PointProxy {
  Facet_handle seed;
  Point_3 center;
  Vector_3 normal;
};

struct CompactMetric {
  CompactMetric(const FacetCenterMap &_center_pmap)
    : center_pmap(_center_pmap) {}

  FT operator()(const Facet_const_handle &f, const PointProxy &px) {
    return FT(std::sqrt(CGAL::to_double(
      CGAL::squared_distance(center_pmap[f], px.center))));
  }

  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
};

struct PointProxyFitting {
  PointProxyFitting(const FacetCenterMap &_center_pmap,
    const FacetAreaMap &_area_pmap,
    const FacetNormalMap &_normal_pmap)
    : center_pmap(_center_pmap),
    area_pmap(_area_pmap),
    normal_pmap(_normal_pmap),
    error_functor(center_pmap) {}

  template<typename FacetIterator>
  PointProxy operator()(const FacetIterator beg, const FacetIterator end) {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = norm + normal_pmap[*fitr] * area_pmap[*fitr];
    }
    norm = norm / std::sqrt(CGAL::to_double(norm.squared_length()));

    // fitting center
    Vector_3 center = CGAL::NULL_VECTOR;
    FT area(0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      center = center + (center_pmap[*fitr] - CGAL::ORIGIN) * area_pmap[*fitr];
      area += area_pmap[*fitr];
    }
    center = center / area;

    // construct proxy
    PointProxy px;
    px.normal = norm;
    px.center = CGAL::ORIGIN + center;

    // update seed
    px.seed = *beg;
    FT err_min = error_functor(*beg, px);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      FT err = error_functor(*fitr, px);
      if (err < err_min) {
        err_min = err;
        px.seed = *fitr;
      }
    }

    return px;
  }

  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
  const FacetNormalMap normal_pmap;
  CompactMetric error_functor;
};

struct ApproxTrait {
  typedef Kernel GeomTraits;
  typedef PointProxy Proxy;
  typedef CompactMetric ErrorMetric;
  typedef PointProxyFitting ProxyFitting;
  typedef CGAL::PlaneFitting<Polyhedron, FacetAreaMap, FacetNormalMap, VertexPointMap> PlaneFitting;

  ApproxTrait(const Polyhedron &_mesh,
    const VertexPointMap &_point_pmap,
    const FacetCenterMap &_center_pmap,
    const FacetAreaMap &_area_pmap,
    const FacetNormalMap &_normal_pmap)
    : mesh(_mesh), point_pmap(_point_pmap), center_pmap(_center_pmap), area_pmap(_area_pmap), normal_pmap(_normal_pmap) {}

  ErrorMetric construct_fit_error_functor() const {
    return ErrorMetric(center_pmap);
  }

  ProxyFitting construct_proxy_fitting_functor() const {
    return ProxyFitting(center_pmap, area_pmap, normal_pmap);
  }

  PlaneFitting construct_plane_fitting_functor() const {
    return PlaneFitting(mesh, area_pmap, normal_pmap, point_pmap);
  }

  const Polyhedron &mesh;
  const VertexPointMap point_pmap;
  const FacetCenterMap center_pmap;
  const FacetAreaMap area_pmap;
  const FacetNormalMap normal_pmap;
};

int main(int argc, char *argv[])
{
  if (argc < 5)
    return 0;

  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // construct facet normal & area map
  std::map<Facet_const_handle, Vector_3> facet_normals;
  std::map<Facet_const_handle, FT> facet_areas;
  std::map<Facet_const_handle, Point_3> facet_centers;
  for(Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    const Halfedge_const_handle he = fitr->halfedge();
    const Point_3 p1 = he->opposite()->vertex()->point();
    const Point_3 p2 = he->vertex()->point();
    const Point_3 p3 = he->next()->vertex()->point();
    Vector_3 normal = CGAL::unit_normal(p1, p2, p3);
    facet_normals.insert(std::pair<Facet_const_handle, Vector_3>(fitr, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p1, p2, p3))));
    facet_areas.insert(std::pair<Facet_const_handle, FT>(fitr, area));
    facet_centers.insert(std::pair<Facet_const_handle, Point_3>(fitr, CGAL::centroid(p1, p2, p3)));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);
  FacetCenterMap center_pmap(facet_centers);
  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<Polyhedron &>(mesh));

  // create a property-map for segment-ids
  typedef std::map<Facet_const_handle, std::size_t> Facet_id_map;
  Facet_id_map internal_facet_id_map;
  for (Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr)
    internal_facet_id_map.insert(std::pair<Facet_const_handle, std::size_t>(fitr, 0));
  boost::associative_property_map<Facet_id_map> proxy_patch_map(internal_facet_id_map);

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  std::vector<Polyhedron::Vertex_handle> anchor_vtx;
  std::vector<std::vector<std::size_t> > bdrs;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;
  CGAL::vsa_mesh_approximation(init, mesh,
    num_proxies,
    num_iterations,
    proxy_patch_map,
    point_pmap,
    tris,
    anchor_pos,
    anchor_vtx,
    bdrs,
    ApproxTrait(mesh, point_pmap, center_pmap, area_pmap, normal_pmap));

  return EXIT_SUCCESS;
}
