#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H

#include <list>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL
{
template<typename TriangleMesh,
  typename Traits = typename TriangleMesh::Traits>
  struct PlaneProxy
{
  typedef Traits GeomTraits;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Plane_3 Plane;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  face_descriptor seed;
  Vector normal;
  Plane fit_plane;
};

template<typename PlaneProxy,
  typename FacetNormalMap,
  typename FacetAreaMap>
  struct L21Metric
{
  L21Metric(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap),
    area_pmap(area_pmap) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;
  typedef typename FacetAreaMap::key_type face_descriptor;

  FT operator()(const face_descriptor &f, const PlaneProxy &px) {
    Vector v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return area_pmap[f] * scalar_product_functor(v, v);
  }

  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

template<typename PlaneProxy,
  typename FacetAreaMap,
  typename VertexPointMap,
  typename TriangleMesh>
  struct L2Metric
{
  L2Metric(const TriangleMesh &_mesh,
    const FacetAreaMap &_area_pmap,
    const VertexPointMap &_point_pmap)
    : mesh(_mesh), area_pmap(_area_pmap), point_pmap(_point_pmap) {}

  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  FT operator()(const face_descriptor &f, const PlaneProxy &px) {
    halfedge_descriptor he = halfedge(f, mesh);
    const Point &p0 = point_pmap[source(he, mesh)];
    const Point &p1 = point_pmap[target(he, mesh)];
    const Point &p2 = point_pmap[target(next(he, mesh), mesh)];
    FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * area_pmap[f] / FT(6);
  }

  const FacetAreaMap area_pmap;
  const VertexPointMap point_pmap;
  const TriangleMesh &mesh;
};

template<typename PlaneProxy,
  typename L21Metric,
  typename FacetNormalMap,
  typename FacetAreaMap>
  struct L21ProxyFitting
{
  L21ProxyFitting(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap),
    area_pmap(area_pmap),
    error_functor(normal_pmap, area_pmap) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  // Fit and construct a proxy
  template<typename FacetIterator>
  PlaneProxy operator()(const FacetIterator beg, const FacetIterator end) {
    CGAL_assertion(beg != end);

    // fitting normal
    FT area(0);
    Vector norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      area += area_pmap[*fitr];
      norm = sum_functor(norm,
        scale_functor(normal_pmap[*fitr], area_pmap[*fitr]));
    }
    // norm = scale_functor(norm, FT(1.0 / CGAL::to_double(area)));
    norm = scale_functor(norm, FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    PlaneProxy px;
    px.normal = norm;

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
    // more update? this is gonna be less and less generic
    // fit is specific to each kind of proxy

    return px;
  }

  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
  L21Metric error_functor;
};

template<typename PlaneProxy,
  typename ErrorMetric,
  typename TriangleMesh,
  typename VertexPointMap,
  typename FacetAreaMap>
  struct PCAPlaneFitting
{
  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  PCAPlaneFitting(const TriangleMesh &_mesh,
    const VertexPointMap &_point_pmap,
    const FacetAreaMap &_area_pmap)
    : mesh(_mesh),
    point_pmap(_point_pmap),
    error_functor(_mesh, _area_pmap, _point_pmap) {
      GeomTraits traits;
      scale_functor = traits.construct_scaled_vector_3_object();
    }

  template<typename FacetIterator>
  PlaneProxy operator()(const FacetIterator beg, const FacetIterator end) {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, mesh);
      const Point_3 &p0 = point_pmap[source(he, mesh)];
      const Point_3 &p1 = point_pmap[target(he, mesh)];
      const Point_3 &p2 = point_pmap[target(next(he, mesh), mesh)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    PlaneProxy px;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      px.fit_plane,
      CGAL::Dimension_tag<2>());
    typename GeomTraits::Vector_3 normal = px.fit_plane.orthogonal_vector();
    normal = scale_functor(normal, FT(1.0 / std::sqrt(CGAL::to_double(normal.squared_length()))));
    px.normal = normal;
    // l2 metric doesn't care about the normal as much as l21 do

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

  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
  Construct_scaled_vector_3 scale_functor;
  ErrorMetric error_functor;
};

// Bundled l21 approximation traits
template<typename PlaneProxy,
  typename L21ErrorMetric,
  typename L21ProxyFitting,
  typename FacetNormalMap,
  typename FacetAreaMap>
  struct L21ApproximationTrait
{
  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef PlaneProxy Proxy;
  typedef L21ErrorMetric ErrorMetric;
  typedef L21ProxyFitting ProxyFitting;

  L21ApproximationTrait(
    const FacetNormalMap &_facet_normal_map,
    const FacetAreaMap &_facet_area_map)
    : normal_pmap(_facet_normal_map),
    area_pmap(_facet_area_map) {}

  // traits function object form
  // construct error functor
  ErrorMetric construct_fit_error_functor() const {
    return ErrorMetric(normal_pmap, area_pmap);
  }

  // construct proxy fitting functor
  ProxyFitting construct_proxy_fitting_functor() const {
    return ProxyFitting(normal_pmap, area_pmap);
  }

private:
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
};

// Bundled l2 approximation traits
template<typename TriangleMesh,
  typename PlaneProxy,
  typename L2ErrorMetric,
  typename PCAPlaneFitting,
  typename VertexPointMap,
  typename FacetAreaMap>
  struct L2ApproximationTrait
{
public:
  typedef typename PlaneProxy::GeomTraits GeomTraits;
  typedef PlaneProxy Proxy;
  typedef L2ErrorMetric ErrorMetric;
  typedef PCAPlaneFitting ProxyFitting;

  L2ApproximationTrait(
    const TriangleMesh &_mesh,
    const VertexPointMap &_point_pmap,
    const FacetAreaMap &_facet_area_map)
    : mesh(_mesh),
    point_pmap(_point_pmap),
    area_pmap(_facet_area_map) {
  }

  // traits function object form
  // construct error functor
  ErrorMetric construct_fit_error_functor() const {
    return ErrorMetric(mesh, area_pmap, point_pmap);
  }

  // construct proxy fitting functor
  ProxyFitting construct_proxy_fitting_functor() const {
    return ProxyFitting(mesh, point_pmap, area_pmap);
  }

private:
  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
  const FacetAreaMap area_pmap;
};


} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
