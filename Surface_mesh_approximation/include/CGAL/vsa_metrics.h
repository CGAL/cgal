#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/graph/graph_traits.hpp>

#include <boost/unordered_map.hpp>
#include <list>

namespace CGAL {
namespace VSA {

/*!
 * \ingroup PkgTSMA
 * @brief Plane proxy class for the Variational Shape Approximation algorithm.
 *
 * Class containing few proxy parameters.
 * Used as default proxy for the `L21Metric` and `L2Metric`
 *
 * \cgalModels `Proxy`
 *
 * @tparam GeomTraits geometric traits
 */
template <typename GeomTraits>
class Plane_proxy
{
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;

public:
  // The proxy normal used in the `L21Metric`.
  Vector_3 normal;
  // The fitting plane of the proxy used in the `L2Metric`.
  Plane_3 fit_plane;
};

/*!
 * \ingroup PkgTSMA
 * @brief L21 metric class for the Variational Shape Approximation algorithm.
 * It is simply a functor that takes a facet and a proxy, returns the L21 error between them.
 *
 * \cgalModels `ErrorMetric`
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 * @tparam with_area_weighing set true to activate area weighing
 */
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  bool with_area_weighing = true,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::VSA::Plane_proxy<GeomTraits> >
class L21_metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21_metric(const TriangleMesh &tm, const VertexPointMap &vpoint_map) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    
    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(tm));
    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(get(fnormal_map, f), scale_functor(px.normal, FT(-1.0)));
    return get(farea_map, f) * scalar_product_functor(v, v);
  }

private:
  Face_normal_map fnormal_map;
  Face_area_map farea_map;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization for vertex point map
template <typename TriangleMesh,
  bool with_area_weighing,
  typename GeomTraits,
  typename PlaneProxy>
class L21_metric<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  with_area_weighing,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21_metric(const TriangleMesh &tm) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(tm));
    
    // construct internal facet normal & area map
    VertexPointMap vpoint_map = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(get(fnormal_map, f), scale_functor(px.normal, FT(-1.0)));
    return get(farea_map, f) * scalar_product_functor(v, v);
  }

private:
  Face_normal_map fnormal_map;
  Face_area_map farea_map;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization without area weighing
template <typename TriangleMesh,
  typename VertexPointMap,
  typename GeomTraits,
  typename PlaneProxy>
class L21_metric<TriangleMesh,
  VertexPointMap,
  false,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21_metric(const TriangleMesh &tm, const VertexPointMap &vpoint_map) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    
    // construct internal facet normal map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(get(fnormal_map,f), scale_functor(px.normal, FT(-1.0)));
    return scalar_product_functor(v, v);
  }

private:
  Face_normal_map fnormal_map;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization for vertex point map without area weighing
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L21_metric<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  false,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21_metric(const TriangleMesh &tm) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    
    // construct internal facet normal map
    VertexPointMap vpoint_map = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(get(fnormal_map,f), scale_functor(px.normal, FT(-1.0)));
    return scalar_product_functor(v, v);
  }

private:
  Face_normal_map fnormal_map;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/*!
 * \ingroup PkgTSMA
 * @brief L21 proxy fitting class for the Variational Shape Approximation algorithm.
 * It is simply a functor that takes a range of facets, fitting the proxy parameters.
 *
 * \cgalModels `ProxyFitting`
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::VSA::Plane_proxy<GeomTraits> >
class L21_proxy_fitting
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor.
  L21_proxy_fitting(const TriangleMesh &tm, const VertexPointMap &vpoint_map) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(tm));
    
    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns the proxy fitted from the facets from beg to end.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = sum_functor(norm,
        scale_functor(get(fnormal_map, *fitr), get(farea_map, *fitr)));
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    Proxy px;
    px.normal = norm;

    return px;
  }

private:
  Face_normal_map fnormal_map;
  Face_area_map farea_map;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L21_proxy_fitting<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor.
  L21_proxy_fitting(const TriangleMesh &tm) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    fnormal_map = CGAL::internal::add_property(Face_normal_tag("VSA-face_normal"),
      const_cast<TriangleMesh &>(tm));
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(tm));

    // construct internal facet normal & area map
    VertexPointMap vpoint_map = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpoint_map[source(he, tm)];
      const Point_3 &p1 = vpoint_map[target(he, tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, tm), tm)];
      put(fnormal_map, f, CGAL::unit_normal(p0, p1, p2));
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns the proxy fitted from the facets from beg to end.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = sum_functor(norm,
                         scale_functor(get(fnormal_map, *fitr), get(farea_map, *fitr)));
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    Proxy px;
    px.normal = norm;

    return px;
  }

private:
  Face_normal_map fnormal_map;
  Face_area_map farea_map;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/*!
 * \ingroup PkgTSMA
 * @brief L2 metric class for the Variational Shape Approximation algorithm.
 * Functor that takes a facet and a proxy, and returns the L2 error between them.
 *
 * \cgalModels `ErrorMetric`
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::VSA::Plane_proxy<GeomTraits> >
class L2_metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L2_metric(const TriangleMesh &tm_, const VertexPointMap &vpoint_map_)
    : tm(&tm_), vpoint_map(vpoint_map_) {
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(*tm));

    BOOST_FOREACH(face_descriptor f, faces(*tm)) {
      const halfedge_descriptor he = halfedge(f, *tm);
      const Point_3 &p0 = vpoint_map[source(he, *tm)];
      const Point_3 &p1 = vpoint_map[target(he, *tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns l2 fitting error from facet f to proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    halfedge_descriptor he = halfedge(f, *tm);
    const Point_3 &p0 = vpoint_map[source(he, *tm)];
    const Point_3 &p1 = vpoint_map[target(he, *tm)];
    const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
    const FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    const FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    const FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    const FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    const FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    const FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * get(farea_map, f) / FT(6.0);
  }

private:
  const TriangleMesh *tm;
  const VertexPointMap vpoint_map;
  Face_area_map farea_map;
};

// specialization
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L2_metric<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L2_metric(const TriangleMesh &tm_)
    : tm(&tm_), vpoint_map(get(boost::vertex_point, const_cast<TriangleMesh &>(tm_))) {
    farea_map = CGAL::internal::add_property(Face_area_tag("VSA-face_area"),
      const_cast<TriangleMesh &>(*tm));

    BOOST_FOREACH(face_descriptor f, faces(*tm)) {
      const halfedge_descriptor he = halfedge(f, *tm);
      const Point_3 &p0 = vpoint_map[source(he, *tm)];
      const Point_3 &p1 = vpoint_map[target(he, *tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
      put(farea_map, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns l2 fitting error from facet f to proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    halfedge_descriptor he = halfedge(f, *tm);
    const Point_3 &p0 = vpoint_map[source(he, *tm)];
    const Point_3 &p1 = vpoint_map[target(he, *tm)];
    const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
    const FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    const FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    const FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    const FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    const FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    const FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * get(farea_map, f) / FT(6);
  }

private:
  const TriangleMesh *tm;
  const VertexPointMap vpoint_map;
  Face_area_map farea_map;
};

/*!
 * \ingroup PkgTSMA
 * @brief L2 proxy fitting class for the Variational Shape Approximation algorithm.
 * It is simply a functor that takes a range of facets, fitting the L2 proxy parameters.
 * It uses the PCA algorithm to fit the proxy parameters.
 *
 * \cgalModels `ProxyFitting`
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::VSA::Plane_proxy<GeomTraits> >
class L2_proxy_fitting
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // type required by `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // construct L2 proxy fitting functor from a triangle mesh and the vertex point map.
  L2_proxy_fitting(const TriangleMesh &tm_, const VertexPointMap &vpoint_map_)
    : tm(&tm_), vpoint_map(vpoint_map_) {}

  // returns proxy fitted from range of facets.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *tm);
      const Point_3 &p0 = vpoint_map[source(he, *tm)];
      const Point_3 &p1 = vpoint_map[target(he, *tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Proxy px;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      px.fit_plane,
      CGAL::Dimension_tag<2>());

    return px;
  }

private:
  const TriangleMesh *tm;
  const VertexPointMap vpoint_map;
};

// specialization.
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L2_proxy_fitting<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

public:
  // type required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // construct L2 proxy fitting functor from a triangle mesh.
  L2_proxy_fitting(const TriangleMesh &tm_)
    : tm(&tm_), vpoint_map(get(boost::vertex_point, const_cast<TriangleMesh &>(tm_))) {}

  // returns the proxy fitted from a range of facets.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *tm);
      const Point_3 &p0 = vpoint_map[source(he, *tm)];
      const Point_3 &p1 = vpoint_map[target(he, *tm)];
      const Point_3 &p2 = vpoint_map[target(next(he, *tm), *tm)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Proxy px;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      px.fit_plane,
      CGAL::Dimension_tag<2>());

    return px;
  }

private:
  const TriangleMesh *tm;
  const VertexPointMap vpoint_map;
};

} // end namespace VSA
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H
