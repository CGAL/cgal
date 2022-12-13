// Copyright (c) 2015-2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Yin Xu, Andreas Fabri, Ilker O. Yaz
//

#ifndef CGAL_WEIGHTS_PMP_DEPRECATED_H
#define CGAL_WEIGHTS_PMP_DEPRECATED_H

#define CGAL_DEPRECATED_HEADER "<CGAL/Weights/internal/pmp_weights_deprecated.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "This part of the package is deprecated since the version 5.4 of CGAL!"
#include <CGAL/Installation/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

// README:
// This header collects all weights which have been in CGAL before unifying them
// into the new package Weights. This header is for information purpose only. It
// will be removed in the next release.

// Includes.
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {
namespace Weights {
namespace deprecated {

/// \cond SKIP_IN_MANUAL

// Returns the cotangent value of the half angle [v0, v1, v2]
// using the formula in -[Meyer02] Discrete Differential-Geometry Operators for- page 19.
// The potential problem with the previous one (Cotangent_value) is that it does not produce symmetric results
// (i.e. for v0, v1, v2 and v2, v1, v0 the returned cot weights can be slightly different).
// This one provides stable results.
template<typename PolygonMesh>
struct Cotangent_value_Meyer_impl
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  template<typename VertexPointMap>
  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2,
                    VertexPointMap ppmap)
  {
    typedef typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type >::Kernel::Vector_3 Vector;

    const Vector a = get(ppmap, v0) - get(ppmap, v1);
    const Vector b = get(ppmap, v2) - get(ppmap, v1);

    const double dot_ab = to_double(a * b);

    // Rewritten for safer fp arithmetic.
    // double dot_aa = a.squared_length();
    // double dot_bb = b.squared_length();
    // double divider = CGAL::sqrt(dot_aa * dot_bb - dot_ab * dot_ab);

    const Vector cross_ab = CGAL::cross_product(a, b);
    const double divider = CGAL::to_double(CGAL::approximate_sqrt(cross_ab * cross_ab));

    if (divider == 0.0 /* || divider != divider */)
    {
      CGAL::collinear(get(ppmap, v0), get(ppmap, v1), get(ppmap, v2)) ?
        CGAL_warning_msg(false, "Infinite Cotangent value with the degenerate triangle!") :
        CGAL_warning_msg(false, "Infinite Cotangent value due to the floating point arithmetic!");

      return dot_ab > 0.0 ? (std::numeric_limits<double>::max)() :
                            -(std::numeric_limits<double>::max)();
    }

    return dot_ab / divider;
  }
};

// Same as above but with a different API.
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Cotangent_value_Meyer
{
protected:
  typedef VertexPointMap Point_property_map;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;

  const PolygonMesh&  pmesh_;
  Point_property_map ppmap_;

public:
  Cotangent_value_Meyer(const PolygonMesh& pmesh_,
                        VertexPointMap vpmap_)
    : pmesh_(pmesh_), ppmap_(vpmap_)
  { }

  const PolygonMesh& pmesh() { return pmesh_; }
  Point_property_map ppmap() { return ppmap_; }

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    return Cotangent_value_Meyer_impl<PolygonMesh>()(v0, v1, v2, ppmap());
  }
};

// Imported from skeletonization.
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Cotangent_value_Meyer_secure
{
  typedef VertexPointMap Point_property_map;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  PolygonMesh& pmesh_;
  Point_property_map ppmap_;

public:
  Cotangent_value_Meyer_secure(const PolygonMesh& pmesh_,
                               VertexPointMap vpmap_)
    : pmesh_(pmesh_), ppmap_(vpmap_)
  { }

  const PolygonMesh& pmesh() { return pmesh_; }
  Point_property_map ppmap() { return ppmap_; }

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    const Vector a = get(ppmap(), v0) - get(ppmap(), v1);
    const Vector b = get(ppmap(), v2) - get(ppmap(), v1);

    const double dot_ab = CGAL::to_double(a * b);
    const double dot_aa = CGAL::to_double(a.squared_length());
    const double dot_bb = CGAL::to_double(b.squared_length());
    const double lb = -0.999, ub = 0.999;
    double cosine = dot_ab / CGAL::sqrt(dot_aa) / CGAL::sqrt(dot_bb);
    cosine = (cosine < lb) ? lb : cosine;
    cosine = (cosine > ub) ? ub : cosine;
    const double sine = CGAL::sqrt(1.0 - cosine * cosine);
    return cosine / sine;
  }
};

// Returns the cotangent value of the half angle [v0, v1, v2] by clamping between
// [1, 89] degrees as suggested by -[Friedel] Unconstrained Spherical Parameterization-.
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_clamped : CotangentValue
{
  Cotangent_value_clamped() { }

public:
  Cotangent_value_clamped(const PolygonMesh& pmesh_,
                          VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    const double cot_1 = 57.289962;
    const double cot_89 = 0.017455;
    const double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_89, (std::min)(value, cot_1));
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_clamped_2 : CotangentValue
{
  Cotangent_value_clamped_2() { }

public:
  Cotangent_value_clamped_2(const PolygonMesh& pmesh_,
                            VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    const double cot_5 = 5.671282;
    const double cot_175 = -cot_5;
    const double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_175, (std::min)(value, cot_5));
  }
};

template<typename PolygonMesh,
         typename CotangentValue = Cotangent_value_Meyer_impl<PolygonMesh> >
struct Cotangent_value_minimum_zero_impl
  : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  template<class VertexPointMap>
  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2,
                    const VertexPointMap ppmap)
  {
    const double value = CotangentValue::operator()(v0, v1, v2, ppmap);
    return (std::max)(0.0, value);
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_minimum_zero : CotangentValue
{
public:
  Cotangent_value_minimum_zero() { }

  Cotangent_value_minimum_zero(const PolygonMesh& pmesh_,
                               VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    const double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(0.0, value);
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
         typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Voronoi_area
  : CotangentValue
{
public:
  Voronoi_area(const PolygonMesh& pmesh_,
               VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::in_edge_iterator in_edge_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  double operator()(vertex_descriptor v0)
  {
    // return 1.0;
    double voronoi_area = 0.0;
    for (const halfedge_descriptor he :
      halfedges_around_target(halfedge(v0, pmesh()), pmesh()))
    {
      if (is_border(he, pmesh()) ) { continue; }

      CGAL_assertion(CGAL::is_triangle_mesh(pmesh()));
      CGAL_assertion(v0 == target(he, pmesh()));
      const vertex_descriptor v1   = source(he, pmesh());
      const vertex_descriptor v_op = target(next(he, pmesh()), pmesh());

      const Point& v0_p   = get(ppmap(), v0);
      const Point& v1_p   = get(ppmap(), v1);
      const Point& v_op_p = get(ppmap(), v_op);

      // (?) Check if there is a better way to predicate triangle is obtuse or not.
      const CGAL::Angle angle0   = CGAL::angle(v1_p, v0_p, v_op_p);
      const CGAL::Angle angle1   = CGAL::angle(v_op_p, v1_p, v0_p);
      const CGAL::Angle angle_op = CGAL::angle(v0_p, v_op_p, v1_p);

      bool obtuse = (angle0   == CGAL::OBTUSE) ||
                    (angle1   == CGAL::OBTUSE) ||
                    (angle_op == CGAL::OBTUSE);

      if (!obtuse)
      {
        const double cot_v1   = CotangentValue::operator()(v_op, v1, v0);
        const double cot_v_op = CotangentValue::operator()(v0, v_op, v1);

        const double term1 = cot_v1   * to_double((v_op_p - v0_p).squared_length());
        const double term2 = cot_v_op * to_double((v1_p   - v0_p).squared_length());
        voronoi_area += (1.0 / 8.0) * (term1 + term2);

      }
      else
      {
        const double area_t = to_double(CGAL::approximate_sqrt(CGAL::squared_area(v0_p, v1_p, v_op_p)));

        if (angle0 == CGAL::OBTUSE)
          voronoi_area += area_t / 2.0;
        else
          voronoi_area += area_t / 4.0;
      }
    }

    CGAL_warning_msg(voronoi_area != 0.0, "Zero Voronoi area!");
    return voronoi_area;
  }
};

// Returns the cotangent value of the half angle [v0, v1, v2] by dividing the triangle area
// as suggested by -[Mullen08] Spectral Conformal Parameterization-.
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
         typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_area_weighted
  : CotangentValue
{
  Cotangent_value_area_weighted() { }

public:
  Cotangent_value_area_weighted(const PolygonMesh& pmesh_,
                                VertexPointMap vpmap_) :
    CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0,
                    vertex_descriptor v1,
                    vertex_descriptor v2)
  {
    return CotangentValue::operator()(v0, v1, v2) /
            CGAL::sqrt(CGAL::squared_area(get(this->ppmap(), v0),
                                          get(this->ppmap(), v1),
                                          get(this->ppmap(), v2)));
  }
};

// Cotangent weight calculator:
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-.
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-.
template<typename PolygonMesh,
         typename CotangentValue = Cotangent_value_minimum_zero_impl<PolygonMesh> >
struct Cotangent_weight_impl
  : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent weight of the specified halfedge_descriptor.
  // Edge orientation is trivial.
  template<class VertexPointMap>
  double operator()(halfedge_descriptor he,
                    PolygonMesh& pmesh,
                    VertexPointMap ppmap)
  {
    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);

    // Only one triangle for border edges.
    if (is_border_edge(he, pmesh))
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh), pmesh);
      vertex_descriptor v2 = source(he_cw, pmesh);
      if (is_border_edge(he_cw, pmesh))
      {
        const halfedge_descriptor he_ccw = prev(opposite(he, pmesh), pmesh);
        v2 = source(he_ccw, pmesh);
      }

      return (CotangentValue::operator()(v0, v2, v1, ppmap) / 2.0);
    }
    else
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh), pmesh);
      const vertex_descriptor v2 = source(he_cw, pmesh);
      const halfedge_descriptor he_ccw = prev(opposite(he, pmesh), pmesh);
      const vertex_descriptor v3 = source(he_ccw, pmesh);

      return (CotangentValue::operator()(v0, v2, v1, ppmap) / 2.0 +
              CotangentValue::operator()(v0, v3, v1, ppmap) / 2.0 );
    }
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type,
         typename CotangentValue = Cotangent_value_minimum_zero<PolygonMesh, VertexPointMap> >
class Cotangent_weight
  : CotangentValue
{
public:
  Cotangent_weight() { }

  Cotangent_weight(const PolygonMesh& pmesh_,
                   VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  Cotangent_weight(const PolygonMesh& pmesh_)
    : CotangentValue(pmesh_, get(CGAL::vertex_point, pmesh_))
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;

  // Returns the cotangent weight of the specified halfedge_descriptor.
  // Edge orientation is trivial.
  double operator()(halfedge_descriptor he)
  {
    const vertex_descriptor v0 = target(he, pmesh());
    const vertex_descriptor v1 = source(he, pmesh());

    // Only one triangle for border edges.
    if (is_border_edge(he, pmesh()))
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      vertex_descriptor v2 = source(he_cw, pmesh());
      if (is_border_edge(he_cw, pmesh()))
      {
        const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
        v2 = source(he_ccw, pmesh());
      }

      return (CotangentValue::operator()(v0, v2, v1) / 2.0);
    }
    else
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      const vertex_descriptor v2 = source(he_cw, pmesh());
      const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
      const vertex_descriptor v3 = source(he_ccw, pmesh());

      return (CotangentValue::operator()(v0, v2, v1) / 2.0 +
              CotangentValue::operator()(v0, v3, v1) / 2.0 );
    }
  }
};

// Single cotangent from -[Chao10] Simple Geometric Model for Elastic Deformation.
template<typename PolygonMesh,
         typename CotangentValue = Cotangent_value_Meyer_impl<PolygonMesh> >
struct Single_cotangent_weight_impl
  : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle).
  template <class VertexPointMap>
  double operator()(halfedge_descriptor he,
                    PolygonMesh& pmesh,
                    const VertexPointMap& ppmap)
  {
    if (is_border(he, pmesh))
      return 0.0;

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const vertex_descriptor v_op = target(next(he, pmesh), pmesh);
    return CotangentValue::operator()(v0, v_op, v1, ppmap);
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type,
         typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Single_cotangent_weight
  : CotangentValue
{
  Single_cotangent_weight() { }

public:
  Single_cotangent_weight(const PolygonMesh& pmesh_,
                          VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle).
  double operator()(halfedge_descriptor he)
  {
    if (is_border(he, pmesh()))
      return 0.0;

    const vertex_descriptor v0 = target(he, pmesh());
    const vertex_descriptor v1 = source(he, pmesh());
    const vertex_descriptor v_op = target(next(he, pmesh()), pmesh());
    return CotangentValue::operator()(v0, v_op, v1);
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type,
         typename CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_weight_with_triangle_area
  : CotangentValue
{
  typedef PolygonMesh PM;
  typedef VertexPointMap VPMap;
  typedef typename boost::property_traits<VPMap>::value_type Point;

  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

  Cotangent_weight_with_triangle_area() { }

public:
  Cotangent_weight_with_triangle_area(const PolygonMesh& pmesh_,
                                      VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  { }

  const PolygonMesh& pmesh() { return CotangentValue::pmesh(); }
  VertexPointMap ppmap() { return CotangentValue::ppmap(); }

  double operator()(halfedge_descriptor he)
  {
    const vertex_descriptor v0 = target(he, pmesh());
    const vertex_descriptor v1 = source(he, pmesh());

    // Only one triangle for border edges.
    if (is_border_edge(he, pmesh()))
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      vertex_descriptor v2 = source(he_cw, pmesh());
      if (is_border_edge(he_cw, pmesh()))
      {
        const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
        v2 = source(he_ccw, pmesh());
      }

      const Point& v0_p = get(ppmap(), v0);
      const Point& v1_p = get(ppmap(), v1);
      const Point& v2_p = get(ppmap(), v2);
      const double area_t = to_double(CGAL::sqrt(CGAL::squared_area(v0_p, v1_p, v2_p)));
      return (CotangentValue::operator()(v0, v2, v1) / area_t);
    }
    else
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      const vertex_descriptor v2 = source(he_cw, pmesh());
      const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
      const vertex_descriptor v3 = source(he_ccw, pmesh());

      const Point& v0_p = get(ppmap(), v0);
      const Point& v1_p = get(ppmap(), v1);
      const Point& v2_p = get(ppmap(), v2);
      const Point& v3_p = get(ppmap(), v3);
      const double area_t1 = to_double(CGAL::sqrt(CGAL::squared_area(v0_p, v1_p, v2_p)));
      const double area_t2 = to_double(CGAL::sqrt(CGAL::squared_area(v0_p, v1_p, v3_p)));

      return (CotangentValue::operator()(v0, v2, v1) / area_t1 +
              CotangentValue::operator()(v0, v3, v1) / area_t2 );
    }

    return 0.0;
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Mean_value_weight
{
  const PolygonMesh& pmesh_;
  VertexPointMap vpmap_;

public:
  Mean_value_weight(const PolygonMesh& pmesh_,
                    VertexPointMap vpmap)
    : pmesh_(pmesh_), vpmap_(vpmap)
  { }

  const PolygonMesh& pmesh() { return pmesh_; }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;

  // Returns the mean-value coordinate of the specified halfedge_descriptor.
  // Returns different value for different edge orientation (which is a normal
  // behavior according to the formula).
  double operator()(halfedge_descriptor he)
  {
    const vertex_descriptor v0 = target(he, pmesh());
    const vertex_descriptor v1 = source(he, pmesh());
    const Vector vec = get(vpmap_, v0) - get(vpmap_, v1);
    const double norm = CGAL::sqrt(vec.squared_length());

    // Only one triangle for border edges.
    if (is_border_edge(he, pmesh()))
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      vertex_descriptor v2 = source(he_cw, pmesh());
      if (is_border_edge(he_cw, pmesh()))
      {
        const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
        v2 = source(he_ccw, pmesh());
      }

      return (half_tan_value_2(v1, v0, v2) / norm);
    }
    else
    {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh()), pmesh());
      const vertex_descriptor v2 = source(he_cw, pmesh());
      const halfedge_descriptor he_ccw = prev(opposite(he, pmesh()), pmesh());
      const vertex_descriptor v3 = source(he_ccw, pmesh());
      return (half_tan_value_2(v1, v0, v2) / norm +
              half_tan_value_2(v1, v0, v3) / norm);
    }
  }

private:
  // Returns the tangent value of the half angle v0_v1_v2 / 2.
  double half_tan_value(vertex_descriptor v0,
                        vertex_descriptor v1,
                        vertex_descriptor v2)
  {
    const Vector vec0 = get(vpmap_, v1) - get(vpmap_, v2);
    const Vector vec1 = get(vpmap_, v2) - get(vpmap_, v0);
    const Vector vec2 = get(vpmap_, v0) - get(vpmap_, v1);
    const double e0_square = vec0.squared_length();
    const double e1_square = vec1.squared_length();
    const double e2_square = vec2.squared_length();
    const double e0 = CGAL::sqrt(e0_square);
    const double e2 = CGAL::sqrt(e2_square);
    double cos_angle = (e0_square + e2_square - e1_square) / 2.0 / e0 / e2;
    cos_angle = (std::max)(-1.0, (std::min)(1.0, cos_angle)); // clamp into [-1, 1]
    const double angle = acos(cos_angle);
    return (tan(angle / 2.0));
  }

  // My deviation built on Meyer_02.
  double half_tan_value_2(vertex_descriptor v0,
                          vertex_descriptor v1,
                          vertex_descriptor v2)
  {
    const Vector a = get(vpmap_, v0) - get(vpmap_, v1);
    const Vector b = get(vpmap_, v2) - get(vpmap_, v1);
    const double dot_ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const double dot_aa = a.squared_length();
    const double dot_bb = b.squared_length();
    const double dot_aa_bb = dot_aa * dot_bb;

    const double cos_rep = dot_ab;
    const double sin_rep = CGAL::sqrt(dot_aa_bb  - dot_ab * dot_ab);
    const double normalizer = CGAL::sqrt(dot_aa_bb); // |a| * |b|

    // The formula from [Floater04] page 4:
    // tan(Q/2) = (1 - cos(Q)) / sin(Q).
    return (normalizer - cos_rep) / sin_rep;
  }
};

template<typename PolygonMesh,
         typename PrimaryWeight = Cotangent_weight<PolygonMesh>,
         typename SecondaryWeight = Mean_value_weight<PolygonMesh> >
class Hybrid_weight
  : public PrimaryWeight, SecondaryWeight
{
  PrimaryWeight primary;
  SecondaryWeight secondary;

  Hybrid_weight() { }

public:
  Hybrid_weight(const PolygonMesh& pmesh_)
    : primary(pmesh_), secondary(pmesh_)
  { }

  const PolygonMesh& pmesh() { return primary.pmesh(); }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  double operator()(halfedge_descriptor he)
  {
    const double weight = primary(he);
    // if (weight < 0.0) { std::cout << "Negative weight!" << std::endl; }
    return (weight >= 0.0) ? weight : secondary(he);
  }
};

// Trivial uniform weights (created for test purposes).
template<class PolygonMesh>
class Uniform_weight
{
public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  double operator()(halfedge_descriptor /* e */) { return 1.0; }
};

template<class PolygonMesh>
class Scale_dependent_weight_fairing
{
  const PolygonMesh& pmesh_;

public:
  Scale_dependent_weight_fairing(const PolygonMesh& pmesh_)
    : pmesh_(pmesh_)
  { }

  const PolygonMesh& pmesh() { return pmesh_; }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;

  double w_i(vertex_descriptor /* v_i */) { return 1.0; }

  double w_ij(halfedge_descriptor he)
  {
    const Vector v = target(he, pmesh())->point() - source(he, pmesh())->point();
    const double divider = CGAL::sqrt(v.squared_length());
    if (divider == 0.0)
    {
      CGAL_warning_msg(false, "Scale dependent weight - zero length edge.");
      return (std::numeric_limits<double>::max)();
    }

    return 1.0 / divider;
  }
};

template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Cotangent_weight_with_voronoi_area_fairing
{
  typedef PolygonMesh PM;
  typedef VertexPointMap VPMap;
  Voronoi_area<PM, VPMap> voronoi_functor;
  Cotangent_weight<PM, VPMap, Cotangent_value_Meyer<PM, VPMap> > cotangent_functor;

public:
  Cotangent_weight_with_voronoi_area_fairing(PM& pmesh_)
    : voronoi_functor(pmesh_, get(CGAL::vertex_point, pmesh_)),
      cotangent_functor(pmesh_, get(CGAL::vertex_point, pmesh_))
  { }

  Cotangent_weight_with_voronoi_area_fairing(PM& pmesh_,
                                             VPMap vpmap_)
    : voronoi_functor(pmesh_, vpmap_),
      cotangent_functor(pmesh_, vpmap_)
  { }

  PM& pmesh() { return voronoi_functor.pmesh(); }

  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

  double w_i(vertex_descriptor v_i)
  {
    return 0.5 / voronoi_functor(v_i);
  }

  double w_ij(halfedge_descriptor he)
  {
    return cotangent_functor(he) * 2.0;
  }
};

// Cotangent_value_Meyer has been changed to the version:
// Cotangent_value_Meyer_secure to avoid imprecisions from
// the issue #4706 - https://github.com/CGAL/cgal/issues/4706.
template<typename PolygonMesh,
         typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Cotangent_weight_with_voronoi_area_fairing_secure
{
  typedef PolygonMesh PM;
  typedef VertexPointMap VPMap;
  Voronoi_area<PM, VPMap> voronoi_functor;
  Cotangent_weight<PM, VPMap, Cotangent_value_Meyer_secure<PM, VPMap> > cotangent_functor;

public:
  Cotangent_weight_with_voronoi_area_fairing_secure(PM& pmesh_)
    : voronoi_functor(pmesh_, get(CGAL::vertex_point, pmesh_)),
      cotangent_functor(pmesh_, get(CGAL::vertex_point, pmesh_))
  { }

  Cotangent_weight_with_voronoi_area_fairing_secure(PM& pmesh_,
                                                    VPMap vpmap_)
    : voronoi_functor(pmesh_, vpmap_),
      cotangent_functor(pmesh_, vpmap_)
  { }

  PM& pmesh() { return voronoi_functor.pmesh(); }

  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

  double w_i(vertex_descriptor v_i)
  {
    return 0.5 / voronoi_functor(v_i);
  }

  double w_ij(halfedge_descriptor he)
  {
    return cotangent_functor(he) * 2.0;
  }
};

template<class PolygonMesh>
class Uniform_weight_fairing
{
public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  Uniform_weight_fairing(const PolygonMesh&) { }

  double w_ij(halfedge_descriptor /* e */) { return 1.0; }
  double w_i(vertex_descriptor /* v_i */) { return 1.0; }
};

/// \endcond

} // namespace deprecated
} // namespace Weights
} // namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_WEIGHTS_PMP_DEPRECATED_H
