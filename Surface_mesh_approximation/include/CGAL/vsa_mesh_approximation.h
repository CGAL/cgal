#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/internal/Surface_mesh_approximation/VSA.h>
#include <CGAL/property_map.h>

namespace CGAL
{
/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation a triangulated mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam SegmentPropertyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam PointPropertyMap a property map containing the input mesh vertex point map,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           `TriangleMesh::Point_3` as value type
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position
 * @tparam AnchorVertexContainer a container of extracted anchor vertex
 * @tparam BoundaryContainer a container of proxy patch boundary
 * @tparam GeomTraits geometric kernel

 * @param init select seed initialization
 * @param mesh a triangle mesh
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 * @param segment_ids facet proxy patch id property map
 * @param ppmap mesh vertex point property map
 * @param tris approximation indexed triangle soup
 * @param pos anchor position container
 * @param vtx anchor vertex container
 * @param bdrs proxy patch boundary container
 * @param traits kernel traits
 */
template<typename TriangleMesh,
  typename SegmentPropertyMap,
  typename PointPropertyMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename AnchorVertexContainer,
  typename BoundaryContainer,
  typename GeomTraits>
  void vsa_mesh_approximation(
    const int init,
    const TriangleMesh &mesh,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations,
    SegmentPropertyMap segment_ids,
    PointPropertyMap ppmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    AnchorVertexContainer &vtx,
    BoundaryContainer &bdrs,
    GeomTraits traits) {
  // CGAL_precondition(is_pure_triangle(mesh));

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

  // construct facet normal & area map
  std::map<face_descriptor, Vector> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(mesh)) {
    const halfedge_descriptor he = halfedge(f, mesh);
    const Point p1 = get(ppmap, source(he, mesh));
    const Point p2 = get(ppmap, target(he, mesh));
    const Point p3 = get(ppmap, target(next(he, mesh), mesh));
    Vector normal = CGAL::unit_normal(p1, p2, p3);
    // normal = scale_functor(normal,
    //   FT(1.0 / std::sqrt(CGAL::to_double(normal.squared_length()))));
    facet_normals.insert(std::pair<face_descriptor, Vector>(f, normal));

    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p1, p2, p3))));
    // FT area(std::sqrt(CGAL::to_double(area_functor(p1, p2, p3))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  typedef CGAL::PlaneProxy<TriangleMesh, GeomTraits> PlaneProxy;
  typedef CGAL::L21Metric<PlaneProxy, GeomTraits, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::ProxyFitting<PlaneProxy, GeomTraits, L21Metric, FacetNormalMap, FacetAreaMap> ProxyFitting;
  typedef CGAL::internal::VSA<
    TriangleMesh,
    PlaneProxy,
    L21Metric,
    ProxyFitting,
    GeomTraits,
    FacetNormalMap,
    FacetAreaMap,
    PointPropertyMap> VSA;

  VSA algorithm(mesh, ppmap, normal_pmap, area_pmap, traits);

  switch (init) {
    case VSA::RandomInit:
      algorithm.partition(number_of_segments, number_of_iterations, segment_ids);
      break;
    case VSA::IncrementalInit:
      algorithm.partition_incre(number_of_segments, number_of_iterations, segment_ids);
      break;
    case VSA::HierarchicalInit:
      algorithm.partition_hierarchical(number_of_segments, number_of_iterations, segment_ids);
      break;
  }

  algorithm.extract_mesh(segment_ids, tris);
  BOOST_FOREACH(const typename VSA::Anchor &a, algorithm.collect_anchors()) {
    vtx.push_back(a.vtx);
    pos.push_back(a.pos);
  }

  bdrs = algorithm.collect_borders(segment_ids);
}
}

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
