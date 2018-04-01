#ifndef CGAL_L2_METRIC_PLANE_PROXY_H
#define CGAL_L2_METRIC_PLANE_PROXY_H

#include <CGAL/license/Surface_mesh_approximation.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <list>

namespace CGAL {
namespace VSA {

/// \ingroup PkgTSMA
/// @brief Approximation L2 metric of plane proxy.
///
/// \cgalModels `ErrorMetricProxy`
///
/// @tparam TriangleMesh a triangle `FaceListGraph`
/// @tparam VertexPointMap a property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///    as key type, GeomTraits::Point_3 as value type
/// @tparam GeomTraits geometric traits
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits
    = typename TriangleMesh::Traits>
class L2_metric_plane_proxy {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  /// \name Types
  /// @{
  typedef typename GeomTraits::Plane_3 Proxy;
  /// @}

  /// \name Constructor
  /// @{
  /*!
   * @brief Constructor
   * @param tm triangle mesh
   * @param vpmap vertex point map
   */
  L2_metric_plane_proxy(const TriangleMesh &tm, const VertexPointMap &vpmap)
    : m_tm(&tm), m_vpmap(vpmap){
    m_famap = CGAL::internal::add_property(
      Face_area_tag("VSA-face_area"), const_cast<TriangleMesh &>(*m_tm));

    BOOST_FOREACH(face_descriptor f, faces(*m_tm)) {
      const halfedge_descriptor he = halfedge(f, *m_tm);
      const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
      const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
      const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
      put(m_famap, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }
  /// @}

  /// \name Operations
  /*!
   * @brief Computes the L21 error from a face to a proxy, 
   * using integral (closed-form) computation.
   * @param tm input triangle mesh
   * @param f face_descriptor of a face
   * @param px proxy
   * @return computed error
   */
  FT compute_error(const TriangleMesh &tm, const face_descriptor &f, const Proxy &px) const {
    halfedge_descriptor he = halfedge(f, *m_tm);
    const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
    const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
    const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
    const FT sq_d0 = CGAL::squared_distance(p0, px);
    const FT sq_d1 = CGAL::squared_distance(p1, px);
    const FT sq_d2 = CGAL::squared_distance(p2, px);
    const FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    const FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    const FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + 
            d0 * d1 + d1 * d2 + d2 * d0) * get(m_famap, f) / FT(6.0);
  }

  /*!
   * @brief Fits a proxy from a range of faces, in the L2 sense, with an 
   * integral (closed-form) formulation. The best-fit plane passes
   * through the center of mass and is defined by the two principal
   * components of the integral covariance matrix.
   * @tparam FaceRange range of face descriptors, model of Range.
   * @param tm input triangle mesh
   * @param faces the range of faces to be fitted
   * @return fitted proxy
   */
  template <typename FaceRange>
  Proxy fit_proxy(const TriangleMesh &tm, const FaceRange &faces) const {
    CGAL_assertion(!faces.empty());

    std::list<Triangle_3> tris;
    BOOST_FOREACH(const face_descriptor &f, faces) {
      const halfedge_descriptor he = halfedge(f, *m_tm);
      const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
      const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
      const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Proxy proxy;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      proxy,
      CGAL::Dimension_tag<2>());

    // TODO: check and flip plane normal?

    return proxy;
  }
  /// @}

private:
  const TriangleMesh *m_tm;
  const VertexPointMap m_vpmap;
  Face_area_map m_famap;
};

} // namespace VSA
} // namespace CGAL

#endif // CGAL_L2_METRIC_PLANE_PROXY_H
