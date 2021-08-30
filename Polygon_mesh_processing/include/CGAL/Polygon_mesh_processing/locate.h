// Copyright (c) 2014, 2017, 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Stephen Kiazyk
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
#define CGAL_POLYGON_MESH_PROCESSING_LOCATE_H

#include <CGAL/license/Polygon_mesh_processing/locate.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Default.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Random.h>
#include <CGAL/use.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/variant.hpp>

#include <array>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

// Everywhere in this file:
// If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
// such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
// between the coordinates in `bc` and the vertices of the face `f` is the following:
// - `w0` corresponds to `source(halfedge(f, tm), tm)`
// - `w1` corresponds to `target(halfedge(f, tm), tm)`
// - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// The Ray must have the same ambient dimension as the property map's value type (aka, the point type)

template <typename Point,
          int dim = CGAL::Ambient_dimension<Point>::value>
struct Ray_type_selector
{
  typedef typename CGAL::Kernel_traits<Point>::type                Kernel;
  typedef typename Kernel::Ray_2                                   type;
};

template <typename Point>
struct Ray_type_selector<Point, 3>
{
  typedef typename CGAL::Kernel_traits<Point>::type                Kernel;
  typedef typename Kernel::Ray_3                                   type;
};

// Just for convenience
template <typename TriangleMesh,
          typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
struct Location_traits
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type            Point;
  typedef typename internal::Ray_type_selector<Point>::type                      Ray;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type            Geom_traits;
  typedef typename Geom_traits::FT                                               FT;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor            face_descriptor;

  typedef std::array<FT, 3>                                                      Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates>                    Face_location;
};

} // end namespace internal

/// \ingroup PMP_locate_grp
///
/// A variant used in the function `get_descriptor_from_location()`.
template <typename TriangleMesh>
using descriptor_variant = boost::variant<typename boost::graph_traits<TriangleMesh>::vertex_descriptor,
                                          typename boost::graph_traits<TriangleMesh>::halfedge_descriptor,
                                          typename boost::graph_traits<TriangleMesh>::face_descriptor>;

/// \ingroup PMP_locate_grp
///
/// A triplet of coordinates describing the barycentric coordinates of a point
/// with respect to the vertices of a triangular face.
///
/// \sa `Face_location`
template <typename FT>
using Barycentric_coordinates = std::array<FT, 3>;

/// \ingroup PMP_locate_grp
///
/// If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
/// such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
/// between the coordinates in `bc` and the vertices of the face `f` is the following:
///   - `w0` corresponds to `source(halfedge(f, tm), tm)`
///   - `w1` corresponds to `target(halfedge(f, tm), tm)`
///   - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
template <typename TriangleMesh, typename FT>
using Face_location = std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                Barycentric_coordinates<FT> >;

// forward declarations
template <typename FT, typename TriangleMesh>
bool is_in_face(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& loc,
                const TriangleMesh& tm);

template <typename FT, typename TriangleMesh>
descriptor_variant<TriangleMesh>
get_descriptor_from_location(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                             std::array<FT, 3> >& loc,
                             const TriangleMesh& tm);

// end of forward declarations

namespace internal {

template <typename FT, typename TriangleMesh, typename OutputIterator>
OutputIterator
incident_faces(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                               std::array<FT, 3> >& location,
               const TriangleMesh& tm,
               OutputIterator out)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  const descriptor_variant<TriangleMesh> dv = get_descriptor_from_location(location, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    const vertex_descriptor vd = *vd_ptr;
    for(face_descriptor fd : faces_around_target(halfedge(vd, tm), tm))
      *out++ = fd;
  }
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
  {
    const halfedge_descriptor hd = *hd_ptr;
    *out++ = face(hd, tm);
    *out++ = face(opposite(hd, tm), tm);
  }
  else
  {
    const face_descriptor fd = boost::get<face_descriptor>(dv);
    *out++ = fd;
  }

  return out;
}

// Snapping coordinates for robustness
template <typename FT>
bool
snap_coordinates_to_border(std::array<FT, 3>& coords,
                           const FT tolerance = std::numeric_limits<FT>::epsilon())
{
#ifdef CGAL_PMP_LOCATE_DEBUG
  std::cout << "Pre-snapping: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
  std::cout << "Sum: " << coords[0] + coords[1] + coords[2] << std::endl;
  std::cout << "tolerance: " << tolerance << std::endl;
#endif

  // To still keep a sum roughly equals to 1, keep in memory the small changes
  FT residue(0);
  bool snapped = false;

  for(int i=0; i<3; ++i)
  {
    if(CGAL::abs(coords[i]) <= tolerance)
    {
      snapped = true;
      residue += coords[i];
      coords[i] = FT(0);
    }
    else if(CGAL::abs(FT(1) - coords[i]) <= tolerance)
    {
      snapped = true;
      residue -= FT(1) - coords[i];
      coords[i] = FT(1);
    }
  }

  // Dump the residue into one of the barycentric values that is neither 0 nor 1
  for(int i=0; i<3; ++i)
  {
    if(coords[i] != FT(0) && coords[i] != FT(1))
    {
      coords[i] += residue;
      break;
    }
  }

#ifdef CGAL_PMP_LOCATE_DEBUG
  std::cout << "Post-snapping: " << coords[0] << " "
                                 << coords[1] << " "
                                 << coords[2] << std::endl;
  std::cout << "Sum: " << coords[0] + coords[1] + coords[2] << std::endl;
#endif

  return snapped;
}

template <typename FT, typename TriangleMesh>
bool
snap_location_to_border(std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                  std::array<FT, 3> >& loc,
                        const TriangleMesh /*tm*/,
                        const FT tolerance = std::numeric_limits<FT>::epsilon())
{
  return snap_coordinates_to_border(loc.second, tolerance);
}

template <typename K, typename P, int = P::Ambient_dimension::value>
struct Barycentric_coordinate_calculator // 2D version
{
  std::array<typename K::FT, 3>
  operator()(const P& ip, const P& iq, const P& ir, const P& iquery, const K& k) const
  {
    typedef typename K::FT                        FT;
    typedef typename K::Vector_2                  Vector_2;

    typename K::Construct_point_2 cp2 = k.construct_point_2_object();
    typename K::Construct_vector_2 cv2 = k.construct_vector_2_object();
    typename K::Compute_scalar_product_2 csp2 = k.compute_scalar_product_2_object();

    const typename K::Point_2& p = cp2(ip);
    const typename K::Point_2& q = cp2(iq);
    const typename K::Point_2& r = cp2(ir);
    const typename K::Point_2& query = cp2(iquery);

    Vector_2 v0 = cv2(p, q);
    Vector_2 v1 = cv2(p, r);
    Vector_2 v2 = cv2(p, query);

    FT d00 = csp2(v0, v0);
    FT d01 = csp2(v0, v1);
    FT d11 = csp2(v1, v1);
    FT d20 = csp2(v2, v0);
    FT d21 = csp2(v2, v1);

    FT denom = d00 * d11 - d01 * d01;
    CGAL_assertion((d00 * d11 - d01 * d01) != FT(0)); // denom != 0

    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;

    return CGAL::make_array(FT(FT(1) - v - w), v, w);
  }
};

template <typename K, typename P>
struct Barycentric_coordinate_calculator<K, P, 3 /*3D specialization*/>
{
  std::array<typename K::FT, 3>
  operator()(const P& ip, const P& iq, const P& ir, const P& iquery, const K& k) const
  {
    typedef typename K::FT                        FT;
    typedef typename K::Vector_3                  Vector_3;

    typename K::Construct_point_3 cp3 = k.construct_point_3_object();
    typename K::Construct_vector_3 cv3 = k.construct_vector_3_object();
    typename K::Compute_scalar_product_3 csp3 = k.compute_scalar_product_3_object();

    const typename K::Point_3& p = cp3(ip);
    const typename K::Point_3& q = cp3(iq);
    const typename K::Point_3& r = cp3(ir);
    const typename K::Point_3& query = cp3(iquery);

    Vector_3 v0 = cv3(p, q);
    Vector_3 v1 = cv3(p, r);
    Vector_3 v2 = cv3(p, query);

    FT d00 = csp3(v0, v0);
    FT d01 = csp3(v0, v1);
    FT d11 = csp3(v1, v1);
    FT d20 = csp3(v2, v0);
    FT d21 = csp3(v2, v1);

    CGAL_assertion((d00 * d11 - d01 * d01) != FT(0)); // denom != 0
    FT denom_inv = FT(1) / (d00 * d11 - d01 * d01);

    FT v = (d11 * d20 - d01 * d21) * denom_inv;
    FT w = (d00 * d21 - d01 * d20) * denom_inv;

    return CGAL::make_array(FT(FT(1) - v - w), v, w);
  }
};

template <typename K, typename P, int = P::Ambient_dimension::value>
struct Barycentric_point_constructor // 2D version
{
  typedef typename K::FT                            FT;

  P operator()(const P& p, const FT wp, const P& q, const FT wq, const P& r, const FT wr,
               const K& /*k*/) const
  {
    FT sum = wp + wq + wr;
    CGAL_assertion(sum != FT(0));

    // In theory, this should be compute_x_2(compute_point_2(...)) and construct_P() at the end...
    FT x = (wp * p.x() + wq * q.x() + wr * r.x()) / sum;
    FT y = (wp * p.y() + wq * q.y() + wr * r.y()) / sum;

    return P(x, y);
  }
};

template <typename K, typename P>
struct Barycentric_point_constructor<K, P, 3> // 3D version
{
  typedef typename K::FT                            FT;

  P operator()(const P& p, const FT wp, const P& q, const FT wq, const P& r, const FT wr,
               const K& /*k*/) const
  {
    FT sum = wp + wq + wr;
    CGAL_assertion(sum != FT(0));
    FT x = (wp * p.x() + wq * q.x() + wr * r.x()) / sum;
    FT y = (wp * p.y() + wq * q.y() + wr * r.y()) / sum;
    FT z = (wp * p.z() + wq * q.z() + wr * r.z()) / sum;

    return P(x, y, z);
  }
};

} // namespace internal

/// \ingroup PMP_locate_grp
///
/// \brief Given a set of three points and a query point, computes the barycentric
///        coordinates of the query point with respect to the first three points.
///
/// \tparam GeomTraits the type of a geometric traits. Must be a model of `Kernel` and be compatible
///                    with the template parameter `Point`.
/// \tparam Point the type of a geometric 2D or 3D point
///
/// \param p,q,r three points with respect to whom the barycentric coordinates of `query` will be computed
/// \param query the query point whose barycentric coordinates will be computed
/// \param gt an instance of the geometric traits
///
/// \pre `p`, `q`, and `r` are not collinear.
/// \pre `query` lies on the plane defined by `p`, `q`, and `r`.
///
template <typename GeomTraits, typename Point>
std::array<typename GeomTraits::FT, 3>
barycentric_coordinates(const Point& p, const Point& q, const Point& r, const Point& query,
                        const GeomTraits& gt)
{
  internal::Barycentric_coordinate_calculator<GeomTraits, Point> calculator;
  return calculator(p, q, r, query, gt);
}

template <typename Point>
std::array<typename CGAL::Kernel_traits<Point>::type::FT, 3>
barycentric_coordinates(const Point& p, const Point& q, const Point& r, const Point& query)
{
  typedef typename CGAL::Kernel_traits<Point>::type                     Kernel;

  return barycentric_coordinates<Kernel, Point>(p, q, r, query, Kernel());
}

/// \name Random Location Generation
/// @{

/// \ingroup PMP_locate_grp
///
/// \brief returns a random point over the halfedge `hd`, as a location.
///
/// \details The random point is chosen on the halfedge, meaning that all
///          its barycentric coordinates are positive. It is constructed by uniformly generating
///          a value `t` between `0` and `1`  and setting the barycentric coordinates to `t`, `1-t`,
///          and `0` for respetively the source and target of `hd`, and the third vertex.
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param hd a halfedge of `tm`
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
random_location_on_halfedge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                            const TriangleMesh& tm,
                            CGAL::Random& rnd = get_default_random())
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  const int h_id = halfedge_index_in_face(hd, tm);
  const FT t(rnd.uniform_real(0., 1.));

  std::array<FT, 3> coordinates;
  coordinates[h_id] = t;
  coordinates[(h_id+1)%3] = FT(1)-t;
  coordinates[(h_id+2)%3] = FT(0);

  return std::make_pair(face(hd, tm), coordinates);
}

/// \ingroup PMP_locate_grp
///
/// \brief returns a random point over the face `fd`, as a location.
///
/// \details The random point is on the face, meaning that all its barycentric coordinates
///          are positive.  It is constructed by uniformly picking a value `u` between `0` and `1`,
///          a value `v` between `1-u`, and setting the barycentric coordinates to `u`, `v`, and
///          `1-u-v` for respectively the source and target of `halfedge(fd, tm)`, and the third point.
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
random_location_on_face(typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                        const TriangleMesh& tm,
                        CGAL::Random& rnd = get_default_random())
{
  CGAL_USE(tm);
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

  // calling 'rnd.uniform real' with double in case FT comes from an EPECK kernel (which doesn't seem to work too well)
  FT u(rnd.uniform_real(0., 1.));
  FT v(rnd.uniform_real(0., CGAL::to_double(FT(1) - u)));

  return std::make_pair(fd, CGAL::make_array(u, v, FT(FT(1) - u - v)));
}

/// \ingroup PMP_locate_grp
///
/// \brief returns a random point over the mesh `tm`.
///
/// \details The returned location is obtained by choosing a random face of the mesh and
///          a random point on that face. The barycentric coordinates of the point in the face
///          are thus all positive. Note that all faces have the same probability to be chosen.
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
/// \sa `random_location_on_face()`
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
random_location_on_mesh(const TriangleMesh& tm,
                        CGAL::Random& rnd = get_default_random())
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  face_descriptor fd = CGAL::internal::random_face_in_mesh(tm, rnd);
  return random_location_on_face<FT>(fd, tm, rnd);
}

/// @}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns a descriptor to the simplex of smallest dimension
///        on which the point corresponding to the location lies.
///
/// \details In other words:
///          - if the point lies on a vertex, this function returns a `boost::graph_traits<TriangleMesh>::%vertex_descriptor` `v`;
///          - if the point lies on a halfedge, this function returns a `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` `hd`
///            (note that in that case, `loc.first == face(hd, tm)` holds).
///          - otherwise, this function returns a `boost::graph_traits<TriangleMesh>::%face_descriptor`
///            `fd` (equal to `loc.first`).
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `loc` describes the barycentric coordinates of a point that lives within the face (boundary included),
///      meaning the barycentric coordinates are all positive.
///
template <typename FT, typename TriangleMesh>
descriptor_variant<TriangleMesh>
#ifdef DOXYGEN_RUNNING // just for convenience because template alias do not allow template deduction
get_descriptor_from_location(const Face_location<TriangleMesh, FT>& loc,
#else
get_descriptor_from_location(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                             std::array<FT, 3> >& loc,
#endif
                             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef Barycentric_coordinates<FT>                                             Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& bar = loc.second;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());
  CGAL_precondition(is_in_face(loc, tm));

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == FT(1)) // coordinate at target(hd, tm)
      return target(hd, tm);
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == FT(0)) // coordinate at target(hd, tm)
      return prev(hd, tm);
    hd = next(hd, tm);
  }

  return fd;
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location in a face, returns the geometric position described
///        by these coordinates, as a point.
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param loc the location from which a point is constructed
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{If such traits class is provided, its type `FT` must be identical
///                     to the template parameter `FT` of this function.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
/// \returns a point whose type is the same as the value type of the vertex point property map
///          provided by the user or via named parameters, or the internal point map of the mesh `tm`.
///
template <typename FT, typename TriangleMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Point
construct_point(const Face_location<TriangleMesh, FT>& loc,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Point
construct_point(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& loc,
#endif
                const TriangleMesh& tm,
                const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type            Geom_traits;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type            Point;
  typedef typename boost::property_traits<VertexPointMap>::reference             Point_reference;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                    get_const_property_map(boost::vertex_point, tm));
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  halfedge_descriptor hd = halfedge(loc.first, tm);
  const Point_reference p0 = get(vpm, source(hd, tm));
  const Point_reference p1 = get(vpm, target(hd, tm));
  const Point_reference p2 = get(vpm, target(next(hd, tm), tm));

  internal::Barycentric_point_constructor<Geom_traits, Point> bp_constructor;
  return bp_constructor(p0, loc.second[0], p1, loc.second[1], p2, loc.second[2], gt);
}

template <typename FT, typename TriangleMesh>
typename property_map_value<TriangleMesh, boost::vertex_point_t>::type
#ifdef DOXYGEN_RUNNING
construct_point(const Face_location<TriangleMesh, FT>& loc,
#else
construct_point(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& loc,
#endif
                const TriangleMesh& tm)
{
  return construct_point(loc, tm, parameters::all_default());
}

/// \name Location Predicates
/// @{

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns whether the location is on the vertex `vd` or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param vd a vertex of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_on_vertex(const Face_location<TriangleMesh, FT>& loc,
#else
is_on_vertex(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                             std::array<FT, 3> >& loc,
#endif
             const typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;

  if(!is_in_face(loc, tm))
    return false;

  const descriptor_variant<TriangleMesh> dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
    return (vd == *vd_ptr);

  return false;
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns whether this location is on the halfedge `hd` or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param hd a halfedge of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_on_halfedge(const Face_location<TriangleMesh, FT>& loc,
#else
is_on_halfedge(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                               std::array<FT, 3> >& loc,
#endif
               const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  if(!is_in_face(loc, tm))
    return false;

  const descriptor_variant<TriangleMesh> dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
    return (*vd_ptr == source(hd, tm) || *vd_ptr == target(hd, tm));
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
    return (*hd_ptr == hd);

  return false;
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a set of barycentric coordinates, returns whether those barycentric
///        coordinates correspond to a point within the face (boundary included),
///        that is, if all the barycentric coordinates are positive.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param bar an array of barycentric coordinates
/// \param tm a triangulated surface mesh
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_in_face(const Barycentric_coordinates<FT>& bar,
#else
is_in_face(const std::array<FT, 3>& bar,
#endif
           const TriangleMesh& tm)
{
  CGAL_USE(tm);
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  for(int i=0; i<3; ++i)
  {
    // "|| bar[i] > 1." is not needed because if everything is positive and the sum is '1',
    // then each coefficient is below '1'.
    if(bar[i] < FT(0))
      return false;
  }

  return true;
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns whether the location is in the face (boundary included) or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_in_face(const Face_location<TriangleMesh, FT>& loc,
#else
is_in_face(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                           std::array<FT, 3> >& loc,
#endif
           const TriangleMesh& tm)
{
  return is_in_face(loc.second, tm);
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns whether the location is on the boundary of the face or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_on_face_border(const Face_location<TriangleMesh, FT>& loc,
#else
is_on_face_border(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                  std::array<FT, 3> >& loc,
#endif
                  const TriangleMesh& tm)
{
  if(!is_in_face(loc, tm))
    return false;

  const Barycentric_coordinates<FT>& bar = loc.second;

  for(int i=0; i<3; ++i)
    if(bar[i] == FT(0))
      return true;

  return false;
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a location, returns whether the location is on the border of the mesh or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc a location with `loc.first` a face of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
///
template <typename FT, typename TriangleMesh>
bool
#ifdef DOXYGEN_RUNNING
is_on_mesh_border(const Face_location<TriangleMesh, FT>& loc,
#else
is_on_mesh_border(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                  std::array<FT, 3> >& loc,
#endif
                  const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates<FT>& bar = loc.second;

  if(!is_in_face(bar, tm))
    return false;

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == FT(1)) // coordinate at target(hd, tm)
      return bool(CGAL::is_border(target(hd, tm), tm));
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == FT(0)) // coordinate at target(hd, tm)
      return CGAL::is_border(edge(prev(hd, tm), tm), tm);
    hd = next(hd, tm);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

/// @}

/// \name Point Location
/// @{

/// \ingroup PMP_locate_grp
///
/// \brief returns the location of the given vertex `vd` as a location,
///        that is an ordered pair specifying a face incident to `vd`
///        and the barycentric coordinates of the vertex `vd` in that face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param vd a vertex of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `vd` is not an isolated vertex
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
locate_vertex(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
              const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef Face_location<TriangleMesh, FT>                                 Face_location;

  halfedge_descriptor hd = halfedge(vd, tm);

  // Find a real face in case 'hd' is a border halfedge
  for(halfedge_descriptor thd : halfedges_around_target(hd, tm))
  {
    if(!is_border(thd, tm))
    {
      hd = thd;
      break;
    }
  }

  CGAL_postcondition(!CGAL::is_border(hd, tm)); // must find a 'real' face incident to 'vd'

  face_descriptor fd = face(hd, tm);

  CGAL_assertion(target(hd, tm) == vd);
  CGAL_assertion(fd != boost::graph_traits<TriangleMesh>::null_face());

  // isolated vertex
  if(fd == boost::graph_traits<TriangleMesh>::null_face())
    return Face_location();

  FT coords[3] = { FT(0), FT(0), FT(0) };
  hd = next(hd, tm); // so that source(hd, tm) == vd and it's simpler to handle 'index_in_face'
  int halfedge_local_index = halfedge_index_in_face(hd, tm);
  coords[halfedge_local_index] = FT(1);

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \ingroup PMP_locate_grp
///
/// \brief returns the location of a given vertex as a location in `fd`,
///        that is an ordered pair composed of `fd` and of the barycentric coordinates
///        of the vertex in `fd`.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param vd a vertex of `tm` and a vertex of the face `fd`
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
///
/// \pre `fd` is not the null face
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
locate_vertex(const typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
              const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
              const TriangleMesh& tm)
{
  CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

  FT coords[3] = { FT(0), FT(0), FT(0) };
  std::size_t vertex_local_index = vertex_index_in_face(vd, fd, tm);
  coords[vertex_local_index] = FT(1);

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a point described by a halfedge `hd` and a scalar `t`
///        as `p = (1 - t) * source(hd, tm) + t * target(hd, tm)`,
///        returns this location along the given edge as a location, that is
///        an ordered pair specifying a face containing the location and the
///        barycentric coordinates of that location in that face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param hd a halfedge of `tm`
/// \param t the parametric distance of the desired point along `hd`
/// \param tm a triangulated surface mesh
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
locate_on_halfedge(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                   const FT t,
                   const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  face_descriptor fd = face(hd, tm);
  int edge_local_index = halfedge_index_in_face(hd, tm);

  const FT one_minus_t(FT(1) - t);
  FT coords[3] = { FT(0), FT(0), FT(0) };
  coords[edge_local_index] = one_minus_t;
  coords[(edge_local_index + 1) % 3] = t;

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \ingroup PMP_locate_grp
///
/// \brief Given a point `query` and a face `fd` of a triangulated surface mesh,
///        returns this location as a location, that is
///        an ordered pair composed of `fd` and of the barycentric coordinates of
///        `query` with respect to the vertices of `fd`.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh must be a model of `FaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param query a point, whose type is equal to the value type of the vertex point property map
///              (either user-provided via named parameters or the internal point map of the mesh `tm`)
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{If such traits class is provided, its type `FT` must be identical
///                     to the template parameter `FT` of this function.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{snapping_tolerance}
///     \cgalParamDescription{a tolerance value used to snap barycentric coordinates}
///     \cgalParamType{double}
///     \cgalParamDefault{`0`}
///     \cgalParamExtra{Depending on the geometric traits used, the computation of the barycentric coordinates
///                     might be an inexact construction, thus leading to sometimes surprising values
///                     (e.g. a triplet `[0.5, 0.5, -1-e17]` for a point at the middle of an edge).
///                     The coordinates will be snapped towards `0` and `1` if the difference is smaller
///                     than the tolerance value, while still ensuring that the total sum of the coordinates is `1`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre `fd` is not the null face
///
/// \returns a face location. The type `FT` is deduced from the geometric traits, either provided by
///          the user via named parameters (with `geom_traits`) or using `CGAL::Kernel_traits`
///          and the point type of the vertex point property map in use.
///
template <typename TriangleMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Face_location<TriangleMesh, FT>
locate_in_face(const Point& query,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Face_location
locate_in_face(const typename internal::Location_traits<TriangleMesh, NamedParameters>::Point& query,
#endif
               const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
               const TriangleMesh& tm,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type            Point;
  typedef typename boost::property_traits<VertexPointMap>::reference             Point_reference;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type            Geom_traits;
  typedef typename Geom_traits::FT                                               FT;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                    get_const_property_map(boost::vertex_point, tm));
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  FT snap_tolerance = choose_parameter(get_parameter(np, internal_np::snapping_tolerance), 0);

  CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

  vertex_descriptor vd0 = source(halfedge(fd, tm), tm);
  vertex_descriptor vd1 = target(halfedge(fd, tm), tm);
  vertex_descriptor vd2 = target(next(halfedge(fd, tm), tm), tm);

  const Point_reference p0 = get(vpm, vd0);
  const Point_reference p1 = get(vpm, vd1);
  const Point_reference p2 = get(vpm, vd2);

  std::array<FT, 3> coords = barycentric_coordinates<Geom_traits, Point>(p0, p1, p2, query, gt);

  if(snap_tolerance != FT(0) && !is_in_face(coords, tm))
  {
    std::cerr << "Warning: point " << query << " is not in the input face" << std::endl;
    std::cerr << "Coordinates: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

    // Try to to snap the coordinates, hoping the problem is just a -1e-17ish epsilon
    // pushing the coordinates over the edge
    internal::snap_coordinates_to_border(coords, snap_tolerance);
  }

  return std::make_pair(fd, coords);
}

#ifndef DOXYGEN_RUNNING // because this is in the middle of a @{ @} doxygen group
template <typename TriangleMesh>
typename internal::Location_traits<TriangleMesh>::Face_location
locate_in_face(const typename internal::Location_traits<TriangleMesh>::Point& query,
               const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
               const TriangleMesh& tm)
{
  return locate_in_face(query, f, tm, parameters::all_default());
}
#endif

/// \ingroup PMP_locate_grp
///
/// \brief Given a location and a second face adjacent to the first, returns the location of the point in the second face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is the triplet of barycentric coordinates `(w0, w1, w2)`, the correspondance
///          between the coordinates in `bc` and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam FT must be a model of `FieldNumberType`
/// \tparam TriangleMesh must be a model of `FaceGraph`
///
/// \param loc the first location, with `loc.first` being a face of `tm`
/// \param fd the second face, adjacent to `loc.first`
/// \param tm the triangle mesh to which `fd` belongs
///
/// \pre `loc` corresponds to a point that lies on a face incident to both `loc.first` and `fd`.
///
template <typename FT, typename TriangleMesh>
Face_location<TriangleMesh, FT>
#ifdef DOXYGEN_RUNNING
locate_in_adjacent_face(const Face_location<TriangleMesh, FT>& loc,
#else
locate_in_adjacent_face(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                        std::array<FT, 3> >& loc,
#endif
                        const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                        const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  CGAL_assertion_code(typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;)

  if(loc.first == fd)
    return loc;

  Face_location<TriangleMesh, FT> loc_in_fd = std::make_pair(fd, CGAL::make_array(FT(0), FT(0), FT(0)));
  const descriptor_variant<TriangleMesh> dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    int index_of_vd = vertex_index_in_face(*vd_ptr, fd, tm);
    loc_in_fd.second[index_of_vd] = FT(1);
    // Note that the barycentric coordinates were initialized to 0,
    // so the second and third coordinates are already set up properly.
  }
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
  {
    // Note that, here, we know that we are _not_ on a vertex
    const halfedge_descriptor hd = *hd_ptr;
    const halfedge_descriptor opp_hd = opposite(hd, tm);
    CGAL_assertion(face(hd, tm) == loc.first);
    CGAL_assertion(face(opp_hd, tm) == fd);
    CGAL_assertion(loc.first != boost::graph_traits<TriangleMesh>::null_face());
    CGAL_assertion(fd != boost::graph_traits<TriangleMesh>::null_face());

    const int index_of_hd = halfedge_index_in_face(hd, tm);
    const int index_of_opp_hd = halfedge_index_in_face(opp_hd, tm);

    // - Coordinates will be non-null at indices `index_of_hd`
    //   and `index_of_hd + 1` in loc.first.
    // - Coordinates will be non-null at indices `index_of_opp_hd`
    //   and `index_of_opp_hd + 1` in f.
    // - The halfedges `hd` and `opp_hd` have opposite directions.
    loc_in_fd.second[index_of_opp_hd] = loc.second[(index_of_hd + 1)%3];
    loc_in_fd.second[(index_of_opp_hd + 1)%3] = loc.second[index_of_hd];
    // note that the barycentric coordinates were initialized at 0,
    // so the third coordinate is already set up properly
  }
  else
  {
    CGAL_assertion_code(const face_descriptor fd2 = boost::get<face_descriptor>(dv);)
    CGAL_assertion(fd2 != boost::graph_traits<TriangleMesh>::null_face());
    CGAL_assertion(fd2 != fd);

    // Calling this function for a location that is (strictly) in a face but
    // asking for the location in a nearby face is meaningless
    CGAL_assertion(false);
  }

  CGAL_postcondition(loc_in_fd.first == fd);
  return loc_in_fd;
}

// not documenting the next two functions as they are too technical
#ifndef DOXYGEN_RUNNING

// Finding a common face to a location and a point
// - the first location must be known
// - the second must be a point in a face incident to get_descriptor_from_location(known_location)
// note: not returning the query location to emphasis that the known location can change too.
template <typename FT, typename TriangleMesh, typename NamedParameters>
bool
locate_in_common_face(std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& known_location,
                      const typename internal::Location_traits<TriangleMesh, NamedParameters>::Point& query,
                      std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& query_location,
                      const TriangleMesh& tm,
                      const NamedParameters& np,
                      const FT tolerance = std::numeric_limits<FT>::epsilon())
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  descriptor_variant<TriangleMesh> dv = get_descriptor_from_location(known_location, tm);

  bool is_query_location_in_face = false;

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    const vertex_descriptor vd = *vd_ptr;
    halfedge_descriptor hd = halfedge(vd, tm);

    for(face_descriptor fd : faces_around_target(hd, tm))
    {
      if(fd == boost::graph_traits<TriangleMesh>::null_face())
        continue;

      // check if 'query' can be found in that face
      query_location = locate_in_face(query, fd, tm, np);
      internal::snap_location_to_border<FT>(query_location, tm, tolerance); // @tmp keep or not ?

      is_query_location_in_face = is_in_face(query_location, tm);

      if(is_query_location_in_face)
        break;
    }
  }
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
  {
    const halfedge_descriptor hd = *hd_ptr;
    face_descriptor fd = face(hd, tm);

    if(fd != boost::graph_traits<TriangleMesh>::null_face())
    {
      query_location = locate_in_face(query, fd, tm, np);
      internal::snap_location_to_border<FT>(query_location, tm, tolerance); // @tmp keep or not ?
      is_query_location_in_face = is_in_face(query_location, tm);
    }

    if(!is_query_location_in_face)
    {
      fd = face(opposite(hd, tm), tm);
      query_location = locate_in_face(query, fd, tm, np);
      is_query_location_in_face = is_in_face(query_location, tm);
    }
  }
  else
  {
    const face_descriptor fd = boost::get<face_descriptor>(dv);

    CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

    query_location = locate_in_face(query, fd, tm, np);
    internal::snap_location_to_border<FT>(query_location, tm, tolerance); // @tmp keep or not ?
    is_query_location_in_face = is_in_face(query_location, tm);
  }

  // if this is not the same face as for 'known_query', change 'known_location'
  if(is_query_location_in_face && query_location.first != known_location.first)
    known_location = locate_in_adjacent_face(known_location, query_location.first, tm);

  return is_query_location_in_face;
}

// Finding a common face to two locations
// - both locations must be known but can change
template <typename FT, typename TriangleMesh>
bool
locate_in_common_face(std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                               std::array<FT, 3> >& first_location,
                      std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                std::array<FT, 3> >& second_location,
                      const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  // Check that we actually have something to do
  if(first_location.first == second_location.first)
    return true;

  bool is_first_location_on_border = is_on_face_border(first_location, tm);
  bool is_second_location_on_border = is_on_face_border(second_location, tm);

  // We have already checked that they have different faces, if neither are on
  // a border, then it's hopeless
  if(!is_first_location_on_border && !is_second_location_on_border)
    return false;

  // Find a common face in the sets of incident faces of each location
  std::set<face_descriptor> first_incident_faces;
  std::set<face_descriptor> second_incident_faces;

  internal::incident_faces(first_location, tm, std::inserter(first_incident_faces, first_incident_faces.begin()));
  internal::incident_faces(second_location, tm, std::inserter(second_incident_faces, second_incident_faces.begin()));

  typename std::set<face_descriptor>::const_iterator fit = first_incident_faces.begin();
  typename std::set<face_descriptor>::const_iterator fend = first_incident_faces.end();
  typename std::set<face_descriptor>::const_iterator sit = second_incident_faces.begin();
  typename std::set<face_descriptor>::const_iterator send = second_incident_faces.end();

  while(fit!=fend && sit!=send)
  {
    if(*fit == boost::graph_traits<TriangleMesh>::null_face())
      ++fit;
    if(*sit == boost::graph_traits<TriangleMesh>::null_face())
      ++sit;

    if(*fit == *sit)
      break;
    else if(*fit < *sit)
      ++fit;
    else
      ++sit;
  }

  if(fit == fend || sit == send) // no common face...
    return false;

  CGAL_assertion(*fit == *sit);
  face_descriptor common_fd = *fit;
  CGAL_assertion(common_fd != boost::graph_traits<TriangleMesh>::null_face());

  if(first_location.first != common_fd)
    first_location = locate_in_adjacent_face(first_location, common_fd, tm);

  if(second_location.first != common_fd)
    second_location = locate_in_adjacent_face(second_location, common_fd, tm);

  CGAL_postcondition(first_location.first == second_location.first);
  return true;
}
#endif // DOXYGEN_RUNNING

/// @}

namespace internal {

template <typename TriangleMesh,
          typename Point,
          int dim = CGAL::Ambient_dimension<Point>::value>
struct Point_to_Point_3 // 2D case
{
  typedef typename GetGeomTraits<TriangleMesh>::type::Point_3    Point_3;

  Point_3 operator()(const Point& p) const { return Point_3(p.x(), p.y(), 0); }
};

template <typename TriangleMesh>
struct Point_to_Point_3<TriangleMesh,
                        typename GetGeomTraits<TriangleMesh>::type::Point_3,
                        3> // 3D case with nothing to do
{
  typedef typename GetGeomTraits<TriangleMesh>::type::Point_3    Point_3;

  const Point_3& operator()(const Point_3& p) const { return p; }
};

template <typename TriangleMesh, typename Point>
struct Point_to_Point_3<TriangleMesh, Point, 3> // Generic 3D case
{
  typedef typename GetGeomTraits<TriangleMesh>::type::Point_3    Point_3;

  Point_3 operator()(const Point& p) const { return Point_3(p.x(), p.y(), p.z()); }
};

template <typename TriangleMesh>
struct Ray_to_Ray_3 // 2D case
{
  typedef typename GetGeomTraits<TriangleMesh>::type                       Geom_traits;
  typedef typename Geom_traits::Ray_2                                      Ray_2;
  typedef typename Geom_traits::Ray_3                                      Ray_3;
  typedef Point_to_Point_3<TriangleMesh, typename Geom_traits::Point_2>    P2_to_P3;

  Ray_3 operator()(const Ray_2& r) const
  {
    P2_to_P3 to_p3;
    return Ray_3(to_p3(r.source()), to_p3(r.second_point()));
  }

  const Ray_3& operator()(const Ray_3& r) const { return r; }
};

// Readable property map that converts the output of a given vertex point map to a 3D point
template <typename TriangleMesh,
          typename VertexPointMap =
            typename property_map_selector<TriangleMesh, boost::vertex_point_t>::const_type>
struct Point_to_Point_3_VPM
{
private:
  typedef VertexPointMap                                                  VPM;
  typedef Point_to_Point_3_VPM<TriangleMesh, VPM>                         Self;

public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::property_traits<VPM>::value_type                Point;
  typedef Point_to_Point_3<TriangleMesh, Point>                           P_to_P3;

  typedef typename CGAL::Kernel_traits<Point>::Kernel                     K;
  typedef typename K::Point_3                                             Point_3;

  // required typedefs
  typedef typename boost::property_traits<VertexPointMap>::key_type       key_type;
  typedef Point_3                                                         value_type;
  typedef value_type                                                      reference;
  typedef boost::readable_property_map_tag                                category;

  // Constructors
  Point_to_Point_3_VPM() : conv_(), vpm_() { } // required for compilation by AABBtraits
  Point_to_Point_3_VPM(const VertexPointMap vpm) : conv_(), vpm_(vpm) { }
  Point_to_Point_3_VPM(const TriangleMesh& mesh)
    : conv_(), vpm_(get_const_property_map(boost::vertex_point, mesh))
  { }

  // Access
  const P_to_P3& converter() const { return conv_; }
  const VertexPointMap& vpm() const { return vpm_; }

  // get function for property map
  inline friend reference get(const Self& pmap, key_type v) {
    return pmap.converter()(get(pmap.vpm(), v));
  }

private:
  // Can't be const nor references due to AABB_traits, so make sure to use property maps!
  P_to_P3 conv_;
  VertexPointMap vpm_;
};

// Two different functions, because the AABB's traits' VPM must match the passed VPM (that is,
// the original VPM wrapped with P_to_P3_VPM if the VPM's value_type was not Point_3)
template <typename TriangleMesh, typename AABBTraits, typename VPM>
void build_AABB_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree,
                     const VPM& wrapped_vpm,
                     typename std::enable_if<
                       std::is_same<
                         typename AABBTraits::Point_3, typename boost::property_traits<VPM>::value_type
                       >::value>::type* = 0)
{
  typename boost::graph_traits<TriangleMesh>::face_iterator ffirst, fbeyond;
  boost::tie(ffirst, fbeyond) = faces(tm);
  outTree.rebuild(ffirst, fbeyond, tm, wrapped_vpm);
  outTree.build();
}

template <typename TriangleMesh, typename AABBTraits, typename VPM>
void build_AABB_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree,
                     const VPM& vpm,
                     typename std::enable_if<
                       !std::is_same<
                         typename AABBTraits::Point_3, typename boost::property_traits<VPM>::value_type
                       >::value>::type* = 0)
{
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VPM>              Wrapped_VPM;
  const Wrapped_VPM wrapped_vpm(vpm);

  return internal::build_AABB_tree(tm, outTree, wrapped_vpm);
}

} // namespace internal

/// \name Nearest Face Location Queries
/// The following functions can be used to find the closest point on a triangle mesh, given either
/// a point or a ray. This closest point is computed using a `CGAL::AABB_tree`. Users intending
/// to call location functions on more than a single point (or ray) should first compute an AABB
/// tree to store it (otherwise, it will be recomputed every time). Note that since the AABB tree
/// class is a 3D structure, it might be required to wrap your point property map to convert your
/// point type to the 3D point type (i.e., your traits' `%Point_3`) if you are working
/// with a 2D triangle structure.
///
/// @{

/// \ingroup PMP_locate_grp
///
/// \brief creates an AABB tree suitable for use with `locate_with_AABB_tree()`.
///
/// \details This function should first be called by users who intend to locate multiple points:
///          in this case, it is better to first build an AABB tree, and use the function
///          `locate_with_AABB_tree()` that takes as parameter an AABB tree, instead of calling `locate()`
///          multiple times, which will build a new AABB tree on every call.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`
/// \tparam Point3VPM must be a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                   as key type and the \cgal 3D point type (your traits' `%Point_3`) as value type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param tm a triangulated surface mesh
/// \param outTree output parameter that stores the computed `AABB_tree`
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{Must be identical to the traits used in the template parameter of the `AABB_traits`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename Point3VPM, typename NamedParameters>
void
build_AABB_tree(const TriangleMesh& tm,
                AABB_tree<
                  AABB_traits<
#ifdef DOXYGEN_RUNNING
                    Geom_traits,
#else
                    typename GetGeomTraits<TriangleMesh, NamedParameters>::type,
#endif
                    CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM> > >& outTree,
                const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type     VertexPointMap;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  const VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                          get_const_property_map(boost::vertex_point, tm));

  return internal::build_AABB_tree(tm, outTree, vpm);
}

#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh, typename AABBTraits>
void build_AABB_tree(const TriangleMesh& tm, AABB_tree<AABBTraits>& outTree)
{
  return build_AABB_tree(tm, outTree, parameters::all_default());
}
#endif

/// \ingroup PMP_locate_grp
///
/// \brief returns the face location nearest to the given point, as a location.
///
/// Note that it is possible for the triangle mesh to have ambiant dimension `2` (e.g. the mesh
/// is a 2D triangulation, or a CGAL::Surface_mesh<CGAL::Point_2<Kernel> >), as long as an appropriate
/// vertex point property map is passed in the AABB tree, which will convert from 2D to 3D.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`
/// \tparam Point3VPM must be a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                   as key type and the \cgal 3D point type (your traits' `%Point_3`) as value type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param p the point to locate on the input triangulated surface mesh
/// \param tree an AABB tree containing the triangular faces of the input surface mesh to perform the point location with
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{Must be identical to the traits used in the template parameter of the `AABB_traits`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{snapping_tolerance}
///     \cgalParamDescription{a tolerance value used to snap barycentric coordinates}
///     \cgalParamType{double}
///     \cgalParamDefault{`0`}
///     \cgalParamExtra{Depending on the geometric traits used, the computation of the barycentric coordinates
///                     might be an inexact construction, thus leading to sometimes surprising values
///                     (e.g. a triplet `[0.5, 0.5, -1-e17]` for a point at the middle of an edge).
///                     The coordinates will be snapped towards `0` and `1` if the difference is smaller
///                     than the tolerance value, while still ensuring that the total sum of the coordinates is `1`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns a face location. The type `FT` is deduced from the geometric traits, either provided by
///          the user via named parameters (with `geom_traits`) or using `CGAL::Kernel_traits`
///          and the point type of the vertex point property map in use.
///
template <typename TriangleMesh, typename Point3VPM, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Face_location<TriangleMesh, FT>
locate_with_AABB_tree(const Point& p,
                      const AABB_tree<AABB_traits<Geom_traits,
                                      AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM> > >& tree,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Face_location
locate_with_AABB_tree(const typename internal::Location_traits<TriangleMesh, NamedParameters>::Point& p,
                      const AABB_tree<AABB_traits<
                                typename GetGeomTraits<TriangleMesh, NamedParameters>::type,
                                CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM> > >& tree,
#endif
                      const TriangleMesh& tm,
                      const NamedParameters& np)
{
  typedef typename internal::Location_traits<TriangleMesh, NamedParameters>::Point         Point;
  typedef internal::Point_to_Point_3<TriangleMesh, Point>                                  P_to_P3;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                      Geom_traits;
  typedef typename CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM>       Primitive;
  typedef typename CGAL::AABB_traits<Geom_traits, Primitive>                               AABB_traits;

  typedef typename Primitive::Point                                                        Point_3;
  CGAL_static_assertion((std::is_same<Point_3, typename P_to_P3::Point_3>::value));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type            VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>                     WrappedVPM;

  const Point_3& p3 = P_to_P3()(p);
  typename AABB_tree<AABB_traits>::Point_and_primitive_id result = tree.closest_point_and_primitive(p3);

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                     Geom_traits;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  // The VPM might return a point of any dimension, but the AABB tree necl1671essarily returns
  // a Point_3. So, wrap the VPM (again) to give a Point_3. Even if it's already wrapped, we're just
  // forwarding a const& anyway.
  const VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                          get_const_property_map(boost::vertex_point, tm));
  const WrappedVPM wrapped_vpm(vpm);

  return locate_in_face(result.first, result.second, tm, CGAL::parameters::vertex_point_map(wrapped_vpm));
}

#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh, typename AABBTraits>
typename internal::Location_traits<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Location_traits<TriangleMesh>::Point& p,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm)
{
  return locate_with_AABB_tree(p, tree, tm, parameters::all_default());
}
#endif

/// \ingroup PMP_locate_grp
///
/// \brief returns the nearest face location to the given point.
///
/// \details Note that this function will build an AABB tree on each call. If you need
///          to call this function more than once, first use `build_AABB_tree()` to create a
///          an AABB tree that you can store and use the function `locate_with_AABB_tree()`.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param p the point to locate on the input triangulated surface mesh
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{snapping_tolerance}
///     \cgalParamDescription{a tolerance value used to snap barycentric coordinates}
///     \cgalParamType{double}
///     \cgalParamDefault{`0`}
///     \cgalParamExtra{Depending on the geometric traits used, the computation of the barycentric coordinates
///                     might be an inexact construction, thus leading to sometimes surprising values
///                     (e.g. a triplet `[0.5, 0.5, -1-e17]` for a point at the middle of an edge).
///                     The coordinates will be snapped towards `0` and `1` if the difference is smaller
///                     than the tolerance value, while still ensuring that the total sum of the coordinates is `1`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Face_location<TriangleMesh, FT>
locate(const Point& p,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Face_location
locate(const typename internal::Location_traits<TriangleMesh, NamedParameters>::Point& p,
#endif
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  // Wrap the input VPM with a one converting to 3D (costs nothing if the input VPM
  // already has value type Kernel::Point_3)
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type          VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>                   WrappedVPM;
  typedef typename internal::Location_traits<TriangleMesh, NamedParameters>::Point       Intrinsic_point;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                    Geom_traits;

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, WrappedVPM>                   AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Geom_traits, AABB_face_graph_primitive>                      AABB_face_graph_traits;

  typedef internal::Point_to_Point_3<TriangleMesh, Intrinsic_point>                      P_to_P3;
  typedef typename AABB_face_graph_traits::Point_3                                       Point_3;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_static_assertion((std::is_same<Point_3, typename P_to_P3::Point_3>::value));

  const VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                          get_const_property_map(boost::vertex_point, tm));
  const WrappedVPM wrapped_vpm(vpm);

  AABB_tree<AABB_face_graph_traits> tree;
  build_AABB_tree(tm, tree, parameters::vertex_point_map(wrapped_vpm));

  const Point_3& p3 = P_to_P3()(p);
  return locate_with_AABB_tree(p3, tree, tm, parameters::vertex_point_map(wrapped_vpm));
}

#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh>
typename internal::Location_traits<TriangleMesh>::Face_location
locate(const typename property_map_value<TriangleMesh, boost::vertex_point_t>::type& p,
       const TriangleMesh& tm)
{
  return locate(p, tm, parameters::all_default());
}
#endif

/// \ingroup PMP_locate_grp
///
/// \brief returns the face location along `ray` nearest to its source point.
///
/// If the ray does not intersect the mesh, a default constructed location is returned.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam Point3VPM must be a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                   as key type and the \cgal 3D point type (your traits' `%Point_3`) as value type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param ray a ray to intersect with the input triangulated surface mesh
/// \param tree an AABB tree containing the triangular faces of the input surface mesh to perform the point location with
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{Must be identical to the traits used in the template parameter of the `AABB_traits`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{snapping_tolerance}
///     \cgalParamDescription{a tolerance value used to snap barycentric coordinates}
///     \cgalParamType{double}
///     \cgalParamDefault{`0`}
///     \cgalParamExtra{Depending on the geometric traits used, the computation of the barycentric coordinates
///                     might be an inexact construction, thus leading to sometimes surprising values
///                     (e.g. a triplet `[0.5, 0.5, -1-e17]` for a point at the middle of an edge).
///                     The coordinates will be snapped towards `0` and `1` if the difference is smaller
///                     than the tolerance value, while still ensuring that the total sum of the coordinates is `1`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre `ray` is an object with the same ambient dimension as the point type (the value type of the vertex point map).
///
template <typename TriangleMesh, typename Point3VPM, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Face_location<TriangleMesh, FT>
locate_with_AABB_tree(const Ray& ray,
                      const AABB_tree<AABB_traits<Geom_traits, AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM> > >& tree,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Face_location
locate_with_AABB_tree(const typename internal::Location_traits<TriangleMesh, NamedParameters>::Ray& ray,
                      const AABB_tree<
                        CGAL::AABB_traits<
                          typename GetGeomTraits<TriangleMesh, NamedParameters>::type,
                          CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM>
                      > >& tree,
#endif
                      const TriangleMesh& tm,
                      const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                        Geom_traits;

  typedef typename Geom_traits::FT                                                           FT;
  typedef typename Geom_traits::Point_3                                                      Point_3;
  typedef typename Geom_traits::Ray_3                                                        Ray_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor                        face_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type              VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>                       WrappedVPM;
  typedef internal::Ray_to_Ray_3<TriangleMesh>                                               R_to_R3;

  typedef typename CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Point3VPM>         Primitive;
  typedef typename CGAL::AABB_traits<Geom_traits, Primitive>                                 AABB_traits;
  typedef AABB_tree<AABB_traits>                                                             AABB_face_graph_tree;
  typedef typename AABB_face_graph_tree::template Intersection_and_primitive_id<Ray_3>::Type Intersection_type;
  typedef boost::optional<Intersection_type>                                                 Ray_intersection;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  // First, transform the ray into a 3D ray if needed
  Ray_3 ray_3 = R_to_R3()(ray);

  std::vector<Ray_intersection> intersections;
  tree.all_intersections(ray_3, std::back_inserter(intersections));

  bool found = false;
  FT nearest_distance = 0;
  Point_3 nearest_point = CGAL::ORIGIN;
  face_descriptor nearest_face;

  for(std::size_t i = 0; i < intersections.size(); ++i)
  {
    if(intersections[i])
    {
      Point_3* intersection_point = boost::get<Point_3>(&(intersections[i]->first));

      if(intersection_point)
      {
        FT distance = CGAL::squared_distance(*intersection_point, ray_3.source());

        if(!found || distance < nearest_distance)
        {
          found = true;
          nearest_point = *intersection_point;
          nearest_distance = distance;
          nearest_face = intersections[i]->second;
        }
      }
    }
  }

  if(found)
  {
    // wrap the VPM to make sure it is producing 3D points
    const VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                            get_const_property_map(boost::vertex_point, tm));
    WrappedVPM wrapped_vpm(vpm);

    return locate_in_face(nearest_point, nearest_face, tm, CGAL::parameters::vertex_point_map(wrapped_vpm));
  }
  else
    return std::make_pair(boost::graph_traits<TriangleMesh>::null_face(),
                          CGAL::make_array(FT(0), FT(0), FT(0)));
}

#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh, typename AABBTraits>
typename internal::Location_traits<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Location_traits<TriangleMesh>::Ray& ray,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm)
{
  return locate_with_AABB_tree(ray, tree, tm, parameters::all_default());
}
#endif

/// \ingroup PMP_locate_grp
///
/// \brief returns the face location along `ray` nearest to its source point.
///
/// If the ray does not intersect the mesh, a default constructed location is returned.
///
/// \details Note that this function will build an AABB tree on each call. If you need
///          to call this function more than once, use `build_AABB_tree()` to cache a
///          copy of the `AABB_tree`, and use the overloads of this function
///          that accept a reference to an AABB tree as input.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param ray a ray to intersect with the input triangulated surface mesh
/// \param tm the input triangulated surface mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{snapping_tolerance}
///     \cgalParamDescription{a tolerance value used to snap barycentric coordinates}
///     \cgalParamType{double}
///     \cgalParamDefault{`0`}
///     \cgalParamExtra{Depending on the geometric traits used, the computation of the barycentric coordinates
///                     might be an inexact construction, thus leading to sometimes surprising values
///                     (e.g. a triplet `[0.5, 0.5, -1-e17]` for a point at the middle of an edge).
///                     The coordinates will be snapped towards `0` and `1` if the difference is smaller
///                     than the tolerance value, while still ensuring that the total sum of the coordinates is `1`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre `ray` is an object with the same ambient dimension as the point type (the value type of the vertex point map).
///
template <typename TriangleMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Face_location<TriangleMesh, FT>
locate(const Ray& ray,
#else
typename internal::Location_traits<TriangleMesh>::Face_location
locate(const typename internal::Location_traits<TriangleMesh, NamedParameters>::Ray& ray,
#endif
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type     VertexPointMap;

  // Wrap the input VPM with a one converting to 3D (costs nothing if the input VPM
  // already has value type Geom_traits::Point_3)
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>              VPM;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type               Geom_traits;

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM>                     AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Geom_traits, AABB_face_graph_primitive>                 AABB_face_graph_traits;
  using parameters::get_parameter;
  using parameters::choose_parameter;

  const VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                          get_const_property_map(boost::vertex_point, tm));
  const VPM wrapped_vpm(vpm);

  AABB_tree<AABB_face_graph_traits> tree;
  build_AABB_tree(tm, tree, parameters::vertex_point_map(wrapped_vpm));

  return locate_with_AABB_tree(ray, tree, tm, np);
}

#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh>
typename internal::Location_traits<TriangleMesh>::Face_location
locate(const typename internal::Ray_type_selector<
               typename internal::Location_traits<TriangleMesh>::Point>::type& ray,
       const TriangleMesh& tm)
{
  return locate(ray, tm, parameters::all_default());
}
#endif

/// @}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
