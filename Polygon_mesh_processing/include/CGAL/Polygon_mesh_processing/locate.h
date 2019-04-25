// Copyright (c) 2014, 2017, 2018 INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Default.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Random.h>

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
// such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
// and the vertices of the face `f` is the following:
// - `w0` corresponds to `source(halfedge(f, tm), tm)`
// - `w1` corresponds to `target(halfedge(f, tm), tm)`
// - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

// to not carry 5-lines-long types everywhere
template <typename TriangleMesh,
          typename NamedParameters = Default>
struct Locate_types
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                  descriptor_variant;

  typedef typename std::conditional<
                     std::is_same<
                       NamedParameters, Default>::value,
                       typename boost::property_map<TriangleMesh,
                                                    boost::vertex_point_t>::const_type,
                       typename GetVertexPointMap<TriangleMesh,
                                                  NamedParameters>::const_type
                     >::type                                               VPM;
  typedef typename boost::property_traits<VPM>::value_type                 Point;

  typedef typename CGAL::Kernel_traits<Point>::type                        Kernel;
  typedef typename Kernel::FT                                              FT;

  typedef std::array<FT, 3>                                                Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates>              Face_location;
};

} // end namespace internal

// forward declarations
template <typename TriangleMesh>
bool is_in_face(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                const TriangleMesh& tm);

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::descriptor_variant
get_descriptor_from_location(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                             const TriangleMesh& tm);

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
               typename internal::Locate_types<TriangleMesh>::FT t,
               const TriangleMesh& tm);
// end of forward declarations

namespace internal {

template<typename TriangleMesh, typename OutputIterator>
OutputIterator
incident_faces(const typename internal::Locate_types<TriangleMesh>::Face_location& location,
               const TriangleMesh& tm,
               OutputIterator out)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                  descriptor_variant;

  descriptor_variant dv = get_descriptor_from_location(location, tm);

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
template<typename TriangleMesh>
bool
snap_coordinates_to_border(typename Locate_types<TriangleMesh>::Barycentric_coordinates& coords,
                           const typename Locate_types<TriangleMesh>::FT tolerance =
                             std::numeric_limits<typename Locate_types<TriangleMesh>::FT>::epsilon())
{
  typedef typename internal::Locate_types<TriangleMesh>::FT            FT;

#ifdef CGAL_PMP_LOCATE_DEBUG
  std::cout << "Pre-snapping: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
  std::cout << "Sum: " << coords[0] + coords[1] + coords[2] << std::endl;
  std::cout << "tolerance: " << tolerance << std::endl;
#endif

  // To still keep a sum roughly equals to 1, keep in memory the small changes
  FT residue = 0.;
  bool snapped = false;

  for(int i=0; i<3; ++i)
  {
    if(CGAL::abs(coords[i]) <= tolerance)
    {
      snapped = true;
      residue += coords[i];
      coords[i] = 0.;
    }
    else if(CGAL::abs(1 - coords[i]) <= tolerance)
    {
      snapped = true;
      residue -= 1. - coords[i];
      coords[i] = 1.;
    }
  }

  // Dump the residue into one of the barycentric values that is neither 0 nor 1
  for(int i=0; i<3; ++i)
  {
    if(coords[i] != 0. && coords[i] != 1.)
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

template<typename TriangleMesh>
bool
snap_location_to_border(typename Locate_types<TriangleMesh>::Face_location& loc,
                        const typename Locate_types<TriangleMesh>::FT tolerance =
                          std::numeric_limits<typename Locate_types<TriangleMesh>::FT>::epsilon())
{
  return snap_coordinates_to_border<TriangleMesh>(loc.second, tolerance);
}

template<typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
random_entity_in_range(InputIterator first, InputIterator beyond,
                       CGAL::Random& rnd = get_default_random())
{
  typedef typename std::iterator_traits<InputIterator>::difference_type  size_type;

  size_type zero = 0, ne = std::distance(first, beyond);
  std::advance(first, rnd.uniform_int(zero, ne - 1));

  return *first;
}

template<typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
random_entity_in_range(const CGAL::Iterator_range<InputIterator>& range,
                       CGAL::Random& rnd = get_default_random())
{
  return random_entity_in_range(range.begin(), range.end(), rnd);
}

template <typename K, typename P, int = P::Ambient_dimension::value>
struct Barycentric_coordinate_calculator // 2D version
{
  CGAL::array<typename K::FT, 3>
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

    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;

    return CGAL::make_array(FT(FT(1) - v - w), v, w);
  }
};

template <typename K, typename P>
struct Barycentric_coordinate_calculator<K, P, 3 /*3D specialization*/>
{
  CGAL::array<typename K::FT, 3>
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

    CGAL_assertion((d00 * d11 - d01 * d01) != FT(0)); // denom != 0.
    FT denom_inv = 1. / (d00 * d11 - d01 * d01);

    FT v = (d11 * d20 - d01 * d21) * denom_inv;
    FT w = (d00 * d21 - d01 * d20) * denom_inv;

    return CGAL::make_array(FT(FT(1) - v - w), v, w);
  }
};

template <typename K, typename P, int = P::Ambient_dimension::value>
struct Barycentric_point_constructor // 2D version
{
  typedef typename K::FT                            FT;

  P operator()(const P& p, const FT wp, const P& q, const FT wq, const P& r, const FT wr) const
  {
    FT sum = wp + wq + wr;
    CGAL_assertion(sum != 0);
    FT x = (wp * p.x() + wq * q.x() + wr * r.x()) / sum;
    FT y = (wp * p.y() + wq * q.y() + wr * r.y()) / sum;

    return P(x, y);
  }
};

template <typename K, typename P>
struct Barycentric_point_constructor<K, P, 3> // 3D version
{
  typedef typename K::FT                            FT;

  P operator()(const P& p, const FT wp, const P& q, const FT wq, const P& r, const FT wr) const
  {
    FT sum = wp + wq + wr;
    CGAL_assertion(sum != 0);
    FT x = (wp * p.x() + wq * q.x() + wr * r.x()) / sum;
    FT y = (wp * p.y() + wq * q.y() + wr * r.y()) / sum;
    FT z = (wp * p.z() + wq * q.z() + wr * r.z()) / sum;

    return P(x, y, z);
  }
};

} // namespace internal

/// \brief returns a random non-null vertex incident to the face `fd` of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `HalfedgeGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
random_vertex_in_face(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                      const PolygonMesh& pm,
                      CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(vertices_around_face(halfedge(fd, pm), pm), rnd);
}

/// \brief returns a random non-null halfedge incident to the face `fd` of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `HalfedgeGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
random_halfedge_in_face(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                        const PolygonMesh& pm,
                        CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(halfedges_around_face(halfedge(fd, pm), pm), rnd);
}

/// \brief returns a random non-null vertex of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `VertexListGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
random_vertex_in_mesh(const PolygonMesh& pm, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(vertices(pm), rnd);
}

/// \brief returns a random non-null halfedge of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `HalfedgeListGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
random_halfedge_in_mesh(const PolygonMesh& pm, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(halfedges(pm), rnd);
}

/// \brief returns a random non-null edge of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `EdgeListGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::edge_descriptor
random_edge_in_mesh(const PolygonMesh& pm, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(edges(pm), rnd);
}

/// \brief returns a random non-null face of the polygon mesh `pm`.
///
/// \tparam PolygonMesh A model of `FaceListGraph`
///
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::face_descriptor
random_face_in_mesh(const PolygonMesh& pm, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(faces(pm), rnd);
}

/// \brief returns the number of `next` one has to apply to the halfedge `hd`
///        for `source(hd, mesh) == vd` to be true, starting from `hd = halfedge(fd, tm)`.
///
/// \tparam PolygonMesh A model of `FaceGraph`
///
/// \param vd a vertex of `pm` whose index we seek
/// \param fd a face of `pm` in which we seek the index of `vd`
/// \param pm a mesh of type `PolygonMesh`
///
/// \pre `vd` is a vertex of `fd`.
template <typename PolygonMesh>
int vertex_index_in_face(const typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                         const typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                         const PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor start = halfedge(fd, pm);
  halfedge_descriptor current = start;
  int counter = 0;

  std::cout << "vd to find: " << &*vd << " pt: " << vd->point() << std::endl;

  do
  {
    std::cout << "current source: " << &*source(current, pm) << " pt: " << source(current, pm)->point() << std::endl;
    if(source(current, pm) == vd)
      break;

    ++counter;
    current = next(current, pm);
  }
  while(current != start);

  CGAL_postcondition(counter == 0 || current != start);

  return counter;
}

/// \brief returns the number of `next` one has to apply to `hd` for `hd == he`
///        to be true, starting from `hd = halfedge(face(hd, tm), tm)`.
///
/// \tparam PolygonMesh A model of `FaceGraph`.
///
/// \param he a halfedge of `pm` whose index in `face(he, tm)` we seek
/// \param pm an object of type `PolygonMesh`
///
template <typename PolygonMesh>
int halfedge_index_in_face(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor he,
                           const PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

  CGAL_precondition(he != boost::graph_traits<PolygonMesh>::null_halfedge());
  CGAL_precondition(!is_border(he, pm));

  face_descriptor f = face(he, pm);
  halfedge_descriptor start = halfedge(f, pm);
  halfedge_descriptor current = start;
  int count = 0;

  while(current != he)
  {
    current = next(current, pm);
    ++count;
  }

  return count;
}

/// \brief Given a set of three points and a query point, computes the barycentric
///        coordinates of the query point with respect to the first three points.
///
/// \pre `p`, `q`, and `r` are not collinear.
/// \pre It must be possible to extract a kernel type model of `Kernel`, using `CGAL::Kernel_traits<P>`
///      (this is the case for all standard %CGAL point types).
/// \pre `query` lies on the plane defined by `p`, `q`, and `r`.
///
template <typename K, typename P>
CGAL::array<typename K::FT, 3>
barycentric_coordinates(const P& p, const P& q, const P& r, const P& query, const K& k)
{
  internal::Barycentric_coordinate_calculator<K, P> calculator;
  return calculator(p, q, r, query, k);
}

template <typename P>
CGAL::array<typename CGAL::Kernel_traits<P>::type::FT, 3>
barycentric_coordinates(const P& p, const P& q, const P& r, const P& query)
{
  typedef typename CGAL::Kernel_traits<P>::type                     Kernel;

  return barycentric_coordinates<Kernel, P>(p, q, r, query, Kernel());
}

// Random locations

/// \brief Returns a random point over the halfedge `hd`, as a `Face_location`.
///
/// \details The point is on the halfedge, meaning that all its barycentric coordinates
///          are positive.
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param hd a halfedge of `tm`
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
random_location_on_halfedge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                            const TriangleMesh& tm,
                            CGAL::Random& rnd = get_default_random())
{
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  FT t(rnd.uniform_real(0., 1.));
  return locate_in_face(hd, t, tm);
}

/// \brief Returns a random point over the face `fd`, as a `Face_location`.
///
/// \details The point is on the face, meaning that all its barycentric coordinates
///          are positive.
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
random_location_on_face(typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                        const TriangleMesh& CGAL_precondition_code(tm),
                        CGAL::Random& rnd = get_default_random())
{
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

  // calling 'rnd.uniform real' with double in case FT comes from an EPECK kernel (which doesn't seem to work too well)
  FT u(rnd.uniform_real(0., 1.));
  FT v(rnd.uniform_real(0., double(FT(1) - u)));

  return std::make_pair(fd, CGAL::make_array(u, v, FT(FT(1) - u - v)));
}

/// \brief Returns a random point over the mesh `tm`.
///
/// \details The returned location corresponds to a random point on a random face
///          of the mesh. The barycentric coordinates of the point in the face
///          are thus all positive.
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param tm a triangulated surface mesh
/// \param rnd optional random number generator
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
random_location_on_mesh(const TriangleMesh& tm, CGAL::Random& rnd = get_default_random())
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  face_descriptor fd = random_face_in_mesh(tm, rnd);
  return random_location_on_face(fd, tm, rnd);
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns a descriptor to the simplex
///        of smallest dimension on which the point corresponding to the location lies.
///
/// \details In other words:
///          - if the point lies on a vertex, this function returns a `vertex_descriptor` `v`;
///          - if the point lies on a halfedge, this function returns a `halfedge_descriptor` `hd`
///            (note that in that case, `loc.first == face(hd, tm)` holds).
///          - otherwise, this function returns a `face_descriptor` `fd` (equal to `loc.first`).
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre the location corresponds to a point that is within a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::descriptor_variant
get_descriptor_from_location(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Barycentric_coordinates  Barycentric_coordinates;

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
    if(bar[i] == 1) // coordinate at target(hd, tm)
      return target(hd, tm);
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == 0) // coordinate at target(hd, tm)
      return prev(hd, tm);
    hd = next(hd, tm);
  }

  return fd;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns the geometric position described
///        by these coordinates, as a point.
///
/// \tparam TriangleMesh A model of `FaceGraph`
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param loc the location to transform into a point
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh, NamedParameters>::Point
location_to_point(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm,
                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type            Point;

  typedef typename internal::Locate_types<TriangleMesh>::Kernel                  Kernel;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_const_property_map(boost::vertex_point, tm));

  halfedge_descriptor hd = halfedge(loc.first, tm);
  const Point& p0 = get(vpm, source(hd, tm));
  const Point& p1 = get(vpm, target(hd, tm));
  const Point& p2 = get(vpm, target(next(hd, tm), tm));

  internal::Barycentric_point_constructor<Kernel, Point> bp_constructor;
  return bp_constructor(p0, loc.second[0], p1, loc.second[1], p2, loc.second[2]);
}

template <typename TriangleMesh>
typename property_map_value<TriangleMesh, boost::vertex_point_t>::type
location_to_point(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  return location_to_point(loc, tm, parameters::all_default());
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns whether the location is
///        on the vertex `vd` or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_on_vertex(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
             const typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::descriptor_variant descriptor_variant;

  if(!is_in_face(loc, tm))
    return false;

  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
    return (v == *vd_ptr);

  return false;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns whether the location is
///        on the halfedge `hd` or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_on_halfedge(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
               const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::descriptor_variant descriptor_variant;

  if(!is_in_face(loc, tm))
    return false;

  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
    return (*vd_ptr == source(h, tm) || *vd_ptr == target(h, tm));
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
    return (*hd_ptr == h);

  return false;
}

/// \brief Given a set of barycentric coordinates, returns whether those barycentric
///        coordinates correspond to a point within the face (boundary included),
///        that is, if all the barycentric coordinates are positive.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_in_face(const typename internal::Locate_types<TriangleMesh>::Barycentric_coordinates& bar,
           const TriangleMesh& CGAL_precondition_code(tm))
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  for(int i=0; i<3; ++i)
  {
    // "|| bar[i] > 1." is not needed because if everything is positive and the sum is '1',
    // then each coefficient is below '1'.
    if(bar[i] < 0.)
      return false;
  }

  return true;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns whether the location is
///        in the face (boundary included) or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_in_face(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
           const TriangleMesh& tm)
{
  return is_in_face(loc.second, tm);
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns whether the location is on the boundary
///        of the face or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_on_face_border(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  typedef typename internal::Locate_types<TriangleMesh>::Face_location    Face_location;
  typedef typename Face_location::second_type                             Barycentric_coordinates;

  if(!is_in_face(loc, tm))
    return false;

  const Barycentric_coordinates& bar = loc.second;

  for(int i=0; i<3; ++i)
    if(bar[i] == 0.)
      return true;

  return false;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, returns whether the location is
///        on the border of the mesh or not.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \pre `loc.first` is a face descriptor corresponding to a face of `tm`.
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
bool
is_on_mesh_border(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Face_location    Face_location;
  typedef typename Face_location::second_type                             Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& bar = loc.second;

  if(!is_in_face(bar, tm))
    return false;

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == 1.) // coordinate at target(hd, tm)
      return bool(CGAL::is_border(target(hd, tm), tm));
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(bar[i] == 0.) // coordinate at target(hd, tm)
      return CGAL::is_border(edge(prev(hd, tm), tm), tm);
    hd = next(hd, tm);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

/// \name Constructions
/// @{

/// \brief Given two faces `fd1` and `fd2` of a polygonal mesh `pm`, returns
///        (if it exists) a common edge. The returned halfedge is incident to `fd1`.
///
/// \tparam TriangleMesh A model of `HalfedgeGraph`
///
template<typename PolygonMesh>
boost::optional<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>
common_halfedge(const typename boost::graph_traits<PolygonMesh>::face_descriptor fd1,
                const typename boost::graph_traits<PolygonMesh>::face_descriptor fd2,
                const PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  if(fd1 == fd2)
    return halfedge(fd1, pm);

  halfedge_descriptor hd = halfedge(fd1, pm), done = hd;
  do
  {
    if(face(opposite(hd, pm), pm) == fd2)
      return hd;

    hd = next(hd, pm);
  }
  while(hd != done);

  return boost::none;
}

/// \name Constructions
/// @{

/// \brief Returns the location of the given vertex `vd` as a `Face_location`,
///        that is an ordered pair specifying a face containing the location
///        and the barycentric coordinates of the vertex `vd` in that face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param vd a vertex of `tm`
/// \param tm a triangulated surface mesh
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  halfedge_descriptor he = halfedge(vd, tm);

  // Find a real face in case 'he' is a border halfedge
  for(halfedge_descriptor hd : halfedges_around_target(he, tm))
  {
    if(!is_border(hd, tm))
    {
      he = hd;
      break;
    }
  }

  CGAL_postcondition(!CGAL::is_border(he, tm)); // must find a 'real' face incident to 'vd'

  face_descriptor fd = face(he, tm);

  CGAL_assertion(target(he, tm) == vd);
  CGAL_assertion(fd != boost::graph_traits<TriangleMesh>::null_face());

  FT coords[3] = { FT(0), FT(0), FT(0) };
  he = next(he, tm); // so that source(he, tm) == vd and it's simpler to handle 'index_in_face'
  std::size_t halfedge_local_index = halfedge_index_in_face(he, tm);
  coords[halfedge_local_index] = FT(1);

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \brief Returns the location of the given vertex as a `Face_location` in `fd`,
///        that is an ordered pair composed of `fd` and of the barycentric coordinates
///        of the vertex in `fd`.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param vd a vertex of `tm` and a vertex of the face `fd`
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(const typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
               const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
               const TriangleMesh& tm)
{
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  FT coords[3] = { FT(0), FT(0), FT(0) };
  std::size_t vertex_local_index = vertex_index_in_face(vd, fd, tm);
  coords[vertex_local_index] = FT(1);

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \brief Given a point described by a halfedge `he` and a scalar `t`
///        as `p = (1 - t) * source(he, tm) + t * target(he, tm)`,
///        returns this location along the given edge as a `Face_location`, that is
///        an ordered pair specifying a face containing the location and the
///        barycentric coordinates of that location in that face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param he a halfedge of `tm`
/// \param t  the parametric distance of the desired point along `he`
/// \param tm a triangulation surface mesh
///
/// \pre `tm` is a triangulated surface mesh.
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
               const typename internal::Locate_types<TriangleMesh>::FT t,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  face_descriptor fd = face(he, tm);
  std::size_t edge_local_index = halfedge_index_in_face(he, tm);

  const FT one_minus_t(FT(1) - t);
  FT coords[3] = { FT(0), FT(0), FT(0) };
  coords[edge_local_index] = one_minus_t;
  coords[(edge_local_index + 1) % 3] = t;

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \brief Given a point `query` and a face `fd` of a triangulation surface mesh,
///        returns this location as a `Face_location`, that is
///        an ordered pair composed of `fd` and of the barycentric coordinates of
///        `query` with respect to the vertices of `fd`.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param query a point
/// \param fd a face of `tm`
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(const typename internal::Locate_types<TriangleMesh, NamedParameters>::Point& query,
               const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
               const TriangleMesh& tm,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Kernel                 K;
  typedef typename K::FT                                                        FT;

  // VertexPointMap
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type           Point;

  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_const_property_map(boost::vertex_point, tm));

  vertex_descriptor vd0 = source(halfedge(fd, tm), tm);
  vertex_descriptor vd1 = target(halfedge(fd, tm), tm);
  vertex_descriptor vd2 = target(next(halfedge(fd, tm), tm), tm);

  const Point& p0 = get(vpm, vd0);
  const Point& p1 = get(vpm, vd1);
  const Point& p2 = get(vpm, vd2);

  std::array<FT, 3> coords = barycentric_coordinates<Point>(p0, p1, p2, query);

  if(!is_in_face(coords, tm))
  {
    std::cerr << "Warning: point " << query << " is not in the input face" << std::endl;
    std::cerr << "Coordinates: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

    // Try to to snap the coordinates, hoping the problem is just a -1e-17ish epsilon
    // pushing the coordinates over the edge
    internal::snap_coordinates_to_border<TriangleMesh>(coords); // @tmp keep or not ?
  }

  return std::make_pair(fd, coords);
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_face(const typename property_map_value<TriangleMesh, boost::vertex_point_t>::type& query,
               const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
               const TriangleMesh& tm)
{
  return locate_in_face(query, f, tm, parameters::all_default());
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        return the `Face_location` of the point in the second face.
///
/// \details If `tm` is the input triangulated surface mesh and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0` corresponds to `source(halfedge(f, tm), tm)`
///          - `w1` corresponds to `target(halfedge(f, tm), tm)`
///          - `w2` corresponds to `target(next(halfedge(f, tm), tm), tm)`
///
/// \tparam TriangleMesh A model of `FaceGraph`
///
/// \param loc the first location, with `loc.first` being a face of `tm`
/// \param fd the second face, adjacent to `loc.first`
/// \param tm the triangle mesh to which `he` belongs
///
/// \pre `loc` corresponds to a point that lies on a face incident to both `loc.first` and `fd`.
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_in_adjacent_face(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                        const typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                        const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  typedef typename internal::Locate_types<TriangleMesh>::Face_location      Face_location;
  typedef typename internal::Locate_types<TriangleMesh>::FT                 FT;

  if(loc.first == fd)
    return loc;

  Face_location loc_in_fd = std::make_pair(fd, CGAL::make_array(FT(0), FT(0), FT(0)));
  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    int index_of_vd = vertex_index_in_face(*vd_ptr, fd, tm);
    loc_in_fd.second[index_of_vd] = 1.;
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
    const face_descriptor fd2 = boost::get<face_descriptor>(dv);
    CGAL_assertion(fd2 != boost::graph_traits<TriangleMesh>::null_face());
    CGAL_assertion(fd2 != fd);

    // Calling this function for a location that is (strictly) in a face but
    // asking for the location in a nearby face is meaningless
    CGAL_assertion(false);
  }

  CGAL_postcondition(loc_in_fd.first == fd);
  return loc_in_fd;
}

// Finding a common face to a location and a point
// - the first location must be known
// - the second must be a point in a face incident to get_descriptor_from_location(known_location)
// note: not returning the query location to emphasis that the known location can change too.
template <typename TriangleMesh>
bool
locate_in_common_face(typename internal::Locate_types<TriangleMesh>::Face_location& known_location,
                      const typename internal::Locate_types<TriangleMesh>::Point& query,
                      typename internal::Locate_types<TriangleMesh>::Face_location& query_location,
                      const TriangleMesh& tm,
                      const typename internal::Locate_types<TriangleMesh>::FT tolerance =
                        std::numeric_limits<typename internal::Locate_types<TriangleMesh>::FT>::epsilon())
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor, halfedge_descriptor, face_descriptor> descriptor_variant;
  descriptor_variant dv = get_descriptor_from_location(known_location, tm);

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
      query_location = locate_in_face(query, fd, tm);
      internal::snap_location_to_border<TriangleMesh>(query_location, tolerance); // @tmp keep or not ?

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
      query_location = locate_in_face(query, fd, tm);
      internal::snap_location_to_border<TriangleMesh>(query_location, tolerance); // @tmp keep or not ?
      is_query_location_in_face = is_in_face(query_location, tm);
    }

    if(!is_query_location_in_face)
    {
      fd = face(opposite(hd, tm), tm);
      query_location = locate_in_face(query, fd, tm);
      is_query_location_in_face = is_in_face(query_location, tm);
    }
  }
  else
  {
    const face_descriptor fd = boost::get<face_descriptor>(dv);

    CGAL_precondition(fd != boost::graph_traits<TriangleMesh>::null_face());

    query_location = locate_in_face(query, fd, tm);
    internal::snap_location_to_border<TriangleMesh>(query_location, tolerance); // @tmp keep or not ?
    is_query_location_in_face = is_in_face(query_location, tm);
  }

  // if this is not the same face as for 'known_query', change 'known_location'
  if(is_query_location_in_face && query_location.first != known_location.first)
    known_location = locate_in_adjacent_face(known_location, query_location.first, tm);

  return is_query_location_in_face;
}

// Finding a common face to two locations
// - both locations must be known but can change
template <typename TriangleMesh>
bool
locate_in_common_face(typename internal::Locate_types<TriangleMesh>::Face_location& first_location,
                      typename internal::Locate_types<TriangleMesh>::Face_location& second_location,
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

/// @}

namespace internal {

template<typename PolygonMesh,
         typename Point,
         int dim = CGAL::Ambient_dimension<Point>::value>
struct Point_to_Point_3 // 2D case
{
  typedef typename internal::Locate_types<PolygonMesh>::Kernel             K;
  typedef typename K::Point_3                                              Point_3;

  Point_3 operator()(const Point& p) const { return Point_3(p.x(), p.y(), 0); }
};

template<typename PolygonMesh>
struct Point_to_Point_3<PolygonMesh,
                        typename internal::Locate_types<PolygonMesh>::Kernel::Point_3,
                        3> // 3D case with nothing to do
{
  typedef typename internal::Locate_types<PolygonMesh>::Kernel             K;
  typedef typename K::Point_3                                              Point_3;

  const Point_3& operator()(const Point_3& p) const { return p; }
};

template<typename PolygonMesh, typename Point>
struct Point_to_Point_3<PolygonMesh, Point, 3> // Generic 3D case
{
  typedef typename internal::Locate_types<PolygonMesh>::Kernel             K;
  typedef typename K::Point_3                                              Point_3;

  Point_3 operator()(const Point& p) const { return Point_3(p.x(), p.y(), p.z()); }
};

template<typename PolygonMesh>
struct Ray_to_Ray_3 // 2D case
{
  typedef typename internal::Locate_types<PolygonMesh>::Kernel             K;
  typedef typename K::Ray_2                                                Ray_2;
  typedef typename K::Ray_3                                                Ray_3;
  typedef Point_to_Point_3<PolygonMesh, typename K::Point_2>               P2_to_P3;

  Ray_3 operator()(const Ray_2& r) const
  {
    P2_to_P3 to_p3;
    return Ray_3(to_p3(r.source()), to_p3(r.second_point()));
  }

  const Ray_3& operator()(const Ray_3& r) const { return r; }
};

/// Readable property map that converts the output of a given vertex point map to a 3D point
template<typename PolygonMesh,
         typename VertexPointMap = typename property_map_selector<PolygonMesh,
                                                                  boost::vertex_point_t>::const_type>
struct Point_to_Point_3_VPM
{
private:
  typedef VertexPointMap                                                  VPM;
  typedef Point_to_Point_3_VPM<PolygonMesh, VPM>                          Self;

public:
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::property_traits<VPM>::value_type                Point;
  typedef Point_to_Point_3<PolygonMesh, Point>                            P_to_P3;

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
  Point_to_Point_3_VPM(const PolygonMesh& mesh)
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
/// @{

/// \brief Creates an `AABB_tree` suitable for use with `locate_with_AABB_tree()`.
///
/// \details This function should be called to create and cache an AABB tree
///          to avoid the AABB tree being rebuilt at every call of `locate_with_AABB_tree()`.
///
/// \tparam TriangleMesh A model of `FaceListGraph`
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param tm a triangulated surface mesh
/// \param outTree output parameter that stores the computed `AABB_tree`
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename AABBTraits, typename NamedParameters>
void build_AABB_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree,
                     const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type     VertexPointMap;

  const VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                                 get_const_property_map(boost::vertex_point, tm));

  return internal::build_AABB_tree(tm, outTree, vpm);
}

template <typename TriangleMesh, typename AABBTraits>
void build_AABB_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree)
{
  return build_AABB_tree(tm, outTree, parameters::all_default());
}

/// \brief Returns the face location nearest to the given point, as a `Face_location`.
///
/// \tparam TriangleMesh A model of `FaceListGraph`
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param p the point to locate on the input triangulated surface mesh
/// \param tree An `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename AABBTraits, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Locate_types<TriangleMesh, NamedParameters>::Point& p,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm,
                      const NamedParameters& np)
{
  typedef typename internal::Locate_types<TriangleMesh, NamedParameters>::Point   Point;
  typedef internal::Point_to_Point_3<TriangleMesh, Point>                         P_to_P3;
  typedef typename AABBTraits::Point_3                                            Point_3;
  CGAL_static_assertion((std::is_same<Point_3, typename P_to_P3::Point_3>::value));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>            WrappedVPM;

  const Point_3& p3 = P_to_P3()(p);
  typename AABB_tree<AABBTraits>::Point_and_primitive_id result = tree.closest_point_and_primitive(p3);

  // The VPM might return a point of any dimension, but the AABB tree necessarily returns
  // a Point_3. So, wrap the VPM (again) to give a Point_3. Even if it's already wrapped, we're just
  // forwarding a const& anyway.
  const VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                                 get_const_property_map(boost::vertex_point, tm));
  const WrappedVPM wrapped_vpm(vpm);

  std::cout << "Result: " << result.first << std::endl;
  std::cout << "in face: " << result.second->vertex(0)->point() << std::endl;
  std::cout << "in face: " << result.second->vertex(1)->point() << std::endl;
  std::cout << "in face: " << result.second->vertex(2)->point() << std::endl;

  return locate_in_face(result.first, result.second, tm, CGAL::parameters::vertex_point_map(wrapped_vpm));
}

template <typename TriangleMesh, typename AABBTraits>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Locate_types<TriangleMesh>::Point& p,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm)
{
  return locate_with_AABB_tree(p, tree, tm, parameters::all_default());
}

/// \brief Returns the nearest face location to the given point.
///
/// \details Note that this function will build an `AABB_tree` on each call. If you need
///          to call this function more than once, use `build_AABB_tree()` to cache a
///          copy of the `AABB_tree`, and use the overloads of this function
///          that accept a reference to an `AABB_tree` as input.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam AABBTraits must be a model of `AABBTraits` used to define a \cgal `AABB_tree`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param p the point to locate on the input triangulated surface mesh
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename internal::Locate_types<TriangleMesh, NamedParameters>::Point& p,
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  // Wrap the input VPM with a one converting to 3D (costs nothing if the input VPM
  // already has value type Kernel::Point_3)
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>            WrappedVPM;
  typedef typename internal::Locate_types<TriangleMesh, NamedParameters>::Point   Intrinsic_point;

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, WrappedVPM>            AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<
            typename internal::Locate_types<TriangleMesh>::Kernel,
            AABB_face_graph_primitive>                                            AABB_face_graph_traits;

  typedef internal::Point_to_Point_3<TriangleMesh, Intrinsic_point>               P_to_P3;
  typedef typename AABB_face_graph_traits::Point_3                                Point_3;

  CGAL_static_assertion((std::is_same<Point_3, typename P_to_P3::Point_3>::value));

  const VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                                 get_const_property_map(boost::vertex_point, tm));
  const WrappedVPM wrapped_vpm(vpm);

  AABB_tree<AABB_face_graph_traits> tree;
  build_AABB_tree(tm, tree, parameters::vertex_point_map(wrapped_vpm));

  const Point_3& p3 = P_to_P3()(p);

  return locate_with_AABB_tree(p3, tree, tm, parameters::vertex_point_map(wrapped_vpm));
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename property_map_value<TriangleMesh, boost::vertex_point_t>::type& p,
       const TriangleMesh& tm)
{
  return locate(p, tm, parameters::all_default());
}

namespace internal {

// The Ray must have the same ambient dimension as the property map's value type (aka, the point type)

template <typename Point,
          int dim = CGAL::Ambient_dimension<Point>::value>
struct Ray_type_selector
{
  typedef typename CGAL::Kernel_traits<Point>::type                Kernel;
  typedef typename Kernel::Ray_2                                   type;
};

template<typename Point>
struct Ray_type_selector<Point, 3>
{
  typedef typename CGAL::Kernel_traits<Point>::type                Kernel;
  typedef typename Kernel::Ray_3                                   type;
};

} // end namespace internal

/// \brief Returns the face location along `ray` nearest to its source point.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam AABBTraits must be a model of `AABBTraits` used to define a \cgal `AABB_tree`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param ray Ray to intersect with the input triangulated surface mesh
/// \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
/// \param tm a triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename AABBTraits, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Ray_type_selector<
                        typename internal::Locate_types<TriangleMesh, NamedParameters>::Point>::type& ray,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm,
                      const NamedParameters& np)
{
  typedef typename internal::Locate_types<TriangleMesh>::Kernel                   Kernel;

  typedef typename Kernel::FT                                                     FT;
  typedef typename Kernel::Point_3                                                Point_3;
  typedef typename Kernel::Ray_3                                                  Ray_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   VertexPointMap;
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>            WrappedVPM;
  typedef internal::Ray_to_Ray_3<TriangleMesh>                                    R_to_R3;

  typedef AABB_tree<AABBTraits>                                                   AABB_face_graph_tree;
  typedef typename AABB_face_graph_tree::template Intersection_and_primitive_id<Ray_3>::Type Intersection_type;
  typedef boost::optional<Intersection_type>                                      Ray_intersection;

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
    const VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                                   get_const_property_map(boost::vertex_point, tm));
    WrappedVPM wrapped_vpm(vpm);

    return locate_in_face(nearest_point, nearest_face, tm, CGAL::parameters::vertex_point_map(wrapped_vpm));
  }
  else
    return std::make_pair(boost::graph_traits<TriangleMesh>::null_face(),
                          CGAL::make_array(FT(0), FT(0), FT(0)));
}

template <typename TriangleMesh, typename AABBTraits>
typename internal::Locate_types<TriangleMesh>::Face_location
locate_with_AABB_tree(const typename internal::Ray_type_selector<
                        typename internal::Locate_types<TriangleMesh>::Point>::type& ray,
                      const AABB_tree<AABBTraits>& tree,
                      const TriangleMesh& tm)
{
  return locate_with_AABB_tree(ray, tree, tm, parameters::all_default());
}

///
/// \brief Returns the face location along `ray` nearest to its source point.
///
/// \details Note that this function will build an `AABB_tree` on each call. If you need
///          to call this function more than once, use `build_AABB_tree()` to cache a
///          copy of the `AABB_tree`, and use the overloads of this function
///          that accept a reference to an `AABB_tree` as input.
///
/// \tparam TriangleMesh must be a model of `FaceListGraph`.
/// \tparam AABBTraits must be a model of `AABBTraits` used to define a \cgal `AABB_tree`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param ray Ray to intersect with the input triangulated surface mesh
/// \param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below:
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map}
///     the property map with the points associated to the vertices of `tm`.
///     If this parameter is omitted, an internal property map for
///     `boost::vertex_point_t` should be available in `TriangleMesh`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename internal::Ray_type_selector<
               typename internal::Locate_types<TriangleMesh, NamedParameters>::Point>::type& ray,
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type     VertexPointMap;

  // Wrap the input VPM with a one converting to 3D (costs nothing if the input VPM
  // already has value type Kernel::Point_3)
  typedef internal::Point_to_Point_3_VPM<TriangleMesh, VertexPointMap>              VPM;

  typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM>                     AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<typename internal::Locate_types<TriangleMesh>::Kernel,
                            AABB_face_graph_primitive>                              AABB_face_graph_traits;

  const VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                                 get_const_property_map(boost::vertex_point, tm));
  const VPM wrapped_vpm(vpm);

  AABB_tree<AABB_face_graph_traits> tree;
  build_AABB_tree(tm, tree, parameters::vertex_point_map(wrapped_vpm));

  return locate_with_AABB_tree(ray, tree, tm, np);
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename internal::Ray_type_selector<
               typename internal::Locate_types<TriangleMesh>::Point>::type& ray,
       const TriangleMesh& tm)
{
  return locate(ray, tm, parameters::all_default());
}

/// @}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
