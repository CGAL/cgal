#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL
{
/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation a triangulated mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 *         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
 *         models of `Hashable`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given
 *         as a named parameter, then the internal one should be initialized
 * @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh with triangulated surface patches to be remeshed
 * @param faces the range of triangular faces defining one or several surface patches to be remeshed
 * @param target_edge_length the edge length that is targetted in the remeshed patch
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of iterations for the
 *    sequence of atomic operations performed (listed in the above description)
 *  \cgalParamEnd
 *  \cgalParamBegin{face_patch_map} a property map with the patch id's associated to the
     faces of `faces`. Instance of a class model of `ReadWritePropertyMap`. It gets
     updated during the remeshing process while new faces are created.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh,
  typename NamedParameters>
bool vsa_mesh_approximation(const TriangleMesh &tm_in,
  TriangleMesh &tm_out,
  const NamedParameters &np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  VPMap point_pmap = choose_param(get_param(np, internal_np::vertex_point),
    get(boost::vertex_point, const_cast<TriangleMesh &>(tm_in)));
    // get_property_map(vertex_point, tm_in));

  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

  typedef CGAL::PlaneProxy<TriangleMesh> PlaneProxy;
  typedef CGAL::L21Metric<TriangleMesh, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::L21ProxyFitting<TriangleMesh, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
  typedef CGAL::VSA_approximation<TriangleMesh, PlaneProxy, L21Metric, L21ProxyFitting> VSAL21;

  // construct facet normal & area map
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(tm_in)) {
    const halfedge_descriptor he = halfedge(f, tm_in);
    const Point_3 p0 = point_pmap[source(he, tm_in)];
    const Point_3 p1 = point_pmap[target(he, tm_in)];
    const Point_3 p2 = point_pmap[target(next(he, tm_in), tm_in)];
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  L21Metric l21_metric(normal_pmap, area_pmap);
  L21ProxyFitting l21_fitting(normal_pmap, area_pmap);

  VSAL21 vsa_l21(l21_metric, l21_fitting);
  vsa_l21.set_mesh(tm_in);

  std::size_t num_proxies = choose_param(get_param(np, internal_np::number_of_proxies),
    num_faces(tm_in) / 100);
  std::size_t num_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 10);
  std::cout << "#px = " << num_proxies << ", #itr = " << num_iterations << std::endl;

  int init = choose_param(get_param(np, internal_np::init_method), 0);
  if (init < 0 || init > 2)
    return false;

  vsa_l21.init_proxies(num_proxies, static_cast<typename VSAL21::Initialization>(init));
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();

  bool pca_plane = choose_param(get_param(np, internal_np::pca_plane) , false);
  FT split_criterion = choose_param(get_param(np, internal_np::chord_subdivide), FT(1));
  bool is_manifold = vsa_l21.meshing(tm_out, split_criterion, pca_plane);

  // vsa_l21.get_proxy_map();

  typedef typename boost::lookup_named_param_def <
      internal_np::anchor_vertex_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<vertex_descriptor> >
    > ::type AnchorVertexOutItr;
  AnchorVertexOutItr avtx_out_itr = choose_param(get_param(np, internal_np::anchor_vertex)
    , std::back_inserter(*(new std::vector<vertex_descriptor>())));
  vsa_l21.get_anchor_vertices(avtx_out_itr);

  typedef typename boost::lookup_named_param_def <
      internal_np::anchor_point_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<Point_3> >
    > ::type AnchorPointOutItr;
  AnchorPointOutItr apts_out_itr = choose_param(get_param(np, internal_np::anchor_point)
    , std::back_inserter(*(new std::vector<Point_3>())));
  vsa_l21.get_anchor_points(apts_out_itr);

  typedef typename boost::lookup_named_param_def <
      internal_np::indexed_triangles_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<int> >
    > ::type IndexedTrisOutItr;
  IndexedTrisOutItr itris_out_itr = choose_param(get_param(np, internal_np::indexed_triangles)
    , std::back_inserter(*(new std::vector<int>())));
  vsa_l21.get_indexed_triangles(itris_out_itr);

  
  // typedef typename boost::lookup_named_param_def <
  //     internal_np::proxies_t,
  //     NamedParameters,
  //     std::back_insert_iterator<std::vector<PlaneProxy> >
  //   > ::type ProxiesOutItr;
  // ProxiesOutItr pxies_out_itr = choose_param(get_param(np, internal_np::proxies)
  //   , std::back_inserter(*(new std::vector<PlaneProxy>())));
  // vsa_l21.get_proxies(pxies_out_itr);

  return is_manifold;
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
