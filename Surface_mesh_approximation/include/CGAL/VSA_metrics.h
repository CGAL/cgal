#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/graph/graph_traits.hpp>

#include <map>
#include <list>

namespace CGAL
{
/*!
 * \ingroup PkgTSMA
 * @brief Plane proxy class for the Variational Shape Approximation algorithm.
 *
 * It is simply a class containing few proxy parameters.
 * It is used as the default proxy for the `L21Metric` and `L2Metric`
 *
 * \cgalModels `Proxy`
 *
 * @tparam GeomTraits geometric traits
 */
template <typename GeomTraits>
class PlaneProxy
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
  typename PlaneProxy = CGAL::PlaneProxy<GeomTraits> >
class L21Metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21Metric(const TriangleMesh &tm, const VertexPointMap &point_pmap)
    : normal_pmap(facet_normals), area_pmap(facet_areas) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return area_pmap[f] * scalar_product_functor(v, v);
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization for vertex point map
template <typename TriangleMesh,
  bool with_area_weighing,
  typename GeomTraits,
  typename PlaneProxy>
class L21Metric<TriangleMesh,
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
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21Metric(const TriangleMesh &tm)
    : normal_pmap(facet_normals), area_pmap(facet_areas) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal & area map
    VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return area_pmap[f] * scalar_product_functor(v, v);
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization without area weighing
template <typename TriangleMesh,
  typename VertexPointMap,
  typename GeomTraits,
  typename PlaneProxy>
class L21Metric<TriangleMesh,
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

  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21Metric(const TriangleMesh &tm, const VertexPointMap &point_pmap)
    : normal_pmap(facet_normals) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return scalar_product_functor(v, v);
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  const FacetNormalMap normal_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization for vertex point map without area weighing
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L21Metric<TriangleMesh,
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
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L21Metric(const TriangleMesh &tm)
    : normal_pmap(facet_normals) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal map
    VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    }
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return scalar_product_functor(v, v);
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  const FacetNormalMap normal_pmap;
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
  typename PlaneProxy = CGAL::PlaneProxy<GeomTraits> >
class L21ProxyFitting
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor.
  L21ProxyFitting(const TriangleMesh &tm, const VertexPointMap &point_pmap)
    : normal_pmap(facet_normals), area_pmap(facet_areas) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
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
        scale_functor(normal_pmap[*fitr], area_pmap[*fitr]));
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    Proxy px;
    px.normal = norm;

    return px;
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

// specialization
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L21ProxyFitting<TriangleMesh,
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
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor.
  L21ProxyFitting(const TriangleMesh &tm)
    : normal_pmap(facet_normals), area_pmap(facet_areas) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();

    // construct internal facet normal & area map
    VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
      facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
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
        scale_functor(normal_pmap[*fitr], area_pmap[*fitr]));
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    Proxy px;
    px.normal = norm;

    return px;
  }

private:
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/*!
 * \ingroup PkgTSMA
 * @brief L2 metric class for the Variational Shape Approximation algorithm.
 * It is simply a functor that takes a facet and a proxy, returns the L2 error between them.
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
  typename PlaneProxy = CGAL::PlaneProxy<GeomTraits> >
class L2Metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L2Metric(const TriangleMesh &tm, const VertexPointMap &_point_pmap)
    : mesh(&tm), area_pmap(facet_areas), point_pmap(_point_pmap) {
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }
  }

  // returns the l2 fitting error of a facet f to proxy px.
  FT operator()(const face_descriptor &f, const PlaneProxy &px) const {
    halfedge_descriptor he = halfedge(f, *mesh);
    const Point_3 &p0 = point_pmap[source(he, *mesh)];
    const Point_3 &p1 = point_pmap[target(he, *mesh)];
    const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
    FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * area_pmap[f] / FT(6);
  }

private:
  const TriangleMesh *mesh;
  std::map<face_descriptor, FT> facet_areas;
  const FacetAreaMap area_pmap;
  const VertexPointMap point_pmap;
};

// specialization
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L2Metric<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // constructor
  L2Metric(const TriangleMesh &tm)
    : mesh(&tm), area_pmap(facet_areas),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(tm))) {
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = point_pmap[source(he, tm)];
      const Point_3 &p1 = point_pmap[target(he, tm)];
      const Point_3 &p2 = point_pmap[target(next(he, tm), tm)];
      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }
  }

  // returns the l2 fitting error of a facet f to proxy px.
  FT operator()(const face_descriptor &f, const PlaneProxy &px) const {
    halfedge_descriptor he = halfedge(f, *mesh);
    const Point_3 &p0 = point_pmap[source(he, *mesh)];
    const Point_3 &p1 = point_pmap[target(he, *mesh)];
    const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
    FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * area_pmap[f] / FT(6);
  }

private:
  const TriangleMesh *mesh;
  std::map<face_descriptor, FT> facet_areas;
  const FacetAreaMap area_pmap;
  const VertexPointMap point_pmap;
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
  typename PlaneProxy = CGAL::PlaneProxy<GeomTraits> >
class L2ProxyFitting
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // construct L2 proxy fitting functor from a triangle mesh and the vertex point map.
  L2ProxyFitting(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap)
    : mesh(&_mesh), point_pmap(_point_pmap) {}

  // returns the proxy fitted from a range of facets.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *mesh);
      const Point_3 &p0 = point_pmap[source(he, *mesh)];
      const Point_3 &p1 = point_pmap[target(he, *mesh)];
      const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
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
  const TriangleMesh *mesh;
  const VertexPointMap point_pmap;
};

// specialization.
template <typename TriangleMesh,
  typename GeomTraits,
  typename PlaneProxy>
class L2ProxyFitting<TriangleMesh,
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  GeomTraits,
  PlaneProxy>
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // construct L2 proxy fitting functor from a triangle mesh.
  L2ProxyFitting(const TriangleMesh &_mesh)
    : mesh(&_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  // returns the proxy fitted from a range of facets.
  template <typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *mesh);
      const Point_3 &p0 = point_pmap[source(he, *mesh)];
      const Point_3 &p1 = point_pmap[target(he, *mesh)];
      const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
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
  const TriangleMesh *mesh;
  const VertexPointMap point_pmap;
};

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_METRICS_H
