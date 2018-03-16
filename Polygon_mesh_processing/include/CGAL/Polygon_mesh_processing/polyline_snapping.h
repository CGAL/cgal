// Copyright (c) 2018 GeometryFactory (France).
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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/tuple.h>

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/value_type.hpp>
#include <boost/foreach.hpp>

#include <cmath>

/// \cond SKIP_IN_MANUAL

namespace CGAL {
namespace Polygon_mesh_processing{
namespace internal {

template <typename Exact_kernel>
bool exact_snapping (const typename Exact_kernel::Segment_3& s0,
                     const typename Exact_kernel::Segment_3& s1,
                     typename Exact_kernel::Point_3& result)
{
  typedef typename Exact_kernel::Point_3 Point_3;
  typedef typename Exact_kernel::Segment_3 Segment_3;
  typedef typename Exact_kernel::Vector_3 Vector_3;
  typedef typename Exact_kernel::Plane_3 Plane_3;

  CGAL_assertion (s0.source() != s0.target());
  CGAL_assertion (s1.source() != s1.target());

  Vector_3 v0 = s0.to_vector();
  Vector_3 v1 = s1.to_vector();
  Vector_3 normal = CGAL::cross_product (v0, v1);
  if (normal == CGAL::NULL_VECTOR)
    return false;
  
  Plane_3 plane0 (s0.source(), normal);
  Plane_3 plane1 (s1.source(), normal);

  Segment_3 s0proj (plane1.projection (s0.source()),
                    plane1.projection (s0.target()));
  Segment_3 s1proj (plane0.projection (s1.source()),
                    plane0.projection (s1.target()));

  CGAL_assertion (s0proj.source() != s0proj.target());
  CGAL_assertion (s1proj.source() != s1proj.target());
  
  if (!CGAL::do_intersect (s0, s1proj) ||
      !CGAL::do_intersect (s1, s0proj))
    return false;

  Point_3 p0, p1;
  
  typename CGAL::cpp11::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection0 = intersection(s0, s1proj);
  if (!intersection0)
    return false;

  if (const Point_3* p = boost::get<Point_3>(&*intersection0))
    p0 = *p;
  else
    return false;

  typename CGAL::cpp11::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection1 = intersection(s1, s0proj);
  if (!intersection1)
    return false;

  if (const Point_3* p = boost::get<Point_3>(&*intersection1))
    p1 = *p;
  else
    return false;

  result = CGAL::midpoint (p0, p1);

  return true;
}

template <class Kernel,
          bool Has_exact_constructions=
          !boost::is_floating_point<typename Kernel::FT>::value >
class Exact_snapping;

template <class Kernel>
class Exact_snapping<Kernel,false>
{
//typedefs
  typedef Kernel Input_kernel;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;
  typedef CGAL::Cartesian_converter<Exact_kernel,Input_kernel> Exact_to_double;

  typedef typename Input_kernel::Point_3 Inexact_point_3;
  typedef typename Input_kernel::Line_3 Inexact_line_3;
  typedef typename Input_kernel::Segment_3 Inexact_segment_3;
  typedef typename Exact_kernel::Point_3 Point_3;
  typedef typename Exact_kernel::Line_3 Line_3;
  typedef typename Exact_kernel::Segment_3 Segment_3;
  
  Exact_kernel ek;
  Exact_to_double exact_to_double;

  Point_3 to_exact(const Inexact_point_3& p) const
  {
    return Point_3(p.x(), p.y(), p.z());
  }
  Inexact_point_3 to_inexact(const Point_3& p, const Exact_to_double& etd) const
  {
    return Inexact_point_3(etd(p.x()),
                           etd(p.y()),
                           etd(p.z()));
  }

public:

  Exact_snapping()
  {

  }

  std::pair<Inexact_point_3, bool> intersection (const Inexact_segment_3& s0, const Inexact_segment_3& s1) const
  {
    Segment_3 exact_s0 (to_exact (s0.source()), to_exact (s0.target()));
    Segment_3 exact_s1 (to_exact (s1.source()), to_exact (s1.target()));
    Point_3 result;
    
    if (exact_snapping<Exact_kernel> (exact_s0, exact_s1, result))
      return std::make_pair (to_inexact (result, exact_to_double), true);

    return std::make_pair (Inexact_point_3(), false);
  }
  
};

template <class Kernel>
class Exact_snapping<Kernel,true>
{
//typedefs
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Segment_3 Segment_3;
public:

  Exact_snapping()
  {

  }

  std::pair<Point_3, bool> intersection (const Segment_3& s0, const Segment_3& s1) const
  {
    Point_3 result;
    bool out = exact_snapping<Kernel> (s0, s1, result);
    return std::make_pair (result, out);
  }
  
};
  
template <typename PolylineRange,
          typename OutputIterator,
          typename Kernel>
struct Intersect_polylines
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::FT FT;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::pair<std::size_t, std::size_t> > Box;
  
  typedef cpp11::tuple<std::size_t, std::size_t, std::size_t, std::size_t, FT, FT> OutputType;
  
  const PolylineRange& polylines;
  double squared_tolerance;
  Exact_snapping<Kernel> exact_snapping;
  mutable OutputIterator output;

  Intersect_polylines (const PolylineRange& polylines,
                       double tolerance,
                       OutputIterator output)
    : polylines (polylines)
    , squared_tolerance (tolerance * tolerance)
    , output (output)
  { }

  void operator()(const Box* b, const Box* c) const
  {
    Segment_3 s1 (polylines[b->info().first][b->info().second],
                  polylines[b->info().first][b->info().second + 1]);
    Segment_3 s2 (polylines[c->info().first][c->info().second],
                  polylines[c->info().first][c->info().second + 1]);
    if (CGAL::squared_distance (s1, s2) > squared_tolerance)
      return;

    Point_3 new_point;
    bool okay;
    boost::tie (new_point, okay) = exact_snapping.intersection (s1, s2);

    if (!okay)
      return;

    FT coord1 = CGAL::approximate_sqrt (CGAL::squared_distance (s1.source(), s1.supporting_line().projection(new_point))
                                        / s1.squared_length());
    FT coord2 = CGAL::approximate_sqrt (CGAL::squared_distance (s2.source(), s2.supporting_line().projection(new_point))
                                        / s2.squared_length());
    
    *(output ++) = cpp11::make_tuple (b->info().first, b->info().second,
                                      c->info().first, c->info().second,
                                      coord1, coord2);
  } 
};


}// namespace internal


// Note this is not officially documented  
/*!
  Computes snapping points among a set of polylines.

  The output is given as a set of `cpp11::tuple` with the following
  elements in the following order:

  - index `i` of the first polyline
  - index `ii` of the segment in polyline `i`
  - index `j` of the second polyline
  - index `jj` of the segment in polyline `j`
  - barycentric coordinate `c_i` of snapping point `p_i` with respect
    to points `p[i][ii]` and `p[i][ii+1]`
  - barycentric coordinate `c_j` of snapping point `p_j` with respect
    to points `p[j][jj]` and `p[j][jj+1]`

  \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange`
  of `Kernel::Point_3`.
  \tparam OutputIterator a model of `OutputIterator` that accepts
  objects of type `cpp11::tuple<std::size_t, std::size_t, std::size_t, std::size_t, FT, FT>`
  \tparam Kernel a model of `Kernel`

  \param polylines the input polyline range.
  \param tolerance the minimum distance between two polylines such
  that the snapping points are computed
  \param output output iterator

 */
template <class PolylineRange, class OutputIterator, class Kernel>
void polyline_snapping (const PolylineRange& polylines,
                        double tolerance,
                        OutputIterator output,
                        const Kernel&)
{
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::pair<std::size_t, std::size_t> > Box;
  typedef typename Kernel::Point_3 Point_3;

  std::size_t size = 0;
  for (std::size_t i = 0; i < polylines.size(); ++ i)
    size += polylines[i].size() - 1;
  
  std::vector<Box> boxes;
  boxes.reserve (size);
  for (std::size_t i = 0; i < polylines.size(); ++ i)
    for (std::size_t j = 0; j < polylines[i].size()-1; ++ j)
    {
      const Point_3& p1 = polylines[i][j];
      const Point_3& p2 = polylines[i][j+1];
      CGAL::Bbox_3 bbox = p1.bbox() + p2.bbox();
      bbox = CGAL::Bbox_3 (bbox.xmin() - tolerance / 2.,
                           bbox.ymin() - tolerance / 2.,
                           bbox.zmin() - tolerance / 2.,
                           bbox.xmax() + tolerance / 2.,
                           bbox.ymax() + tolerance / 2.,
                           bbox.zmax() + tolerance / 2.);
                                 
      boxes.push_back(Box(bbox, std::make_pair (i, j)));
    }

  internal::Intersect_polylines<PolylineRange, OutputIterator, Kernel>
    intersect_polylines(polylines, tolerance, output);
  std::ptrdiff_t cutoff = 2000;
      
  std::size_t idx_start = 0;
  for (std::size_t i = 0; i < polylines.size() - 1; ++ i)
  {
    std::vector<const Box*> box1_ptr
      (boost::make_counting_iterator<const Box*>(&boxes[idx_start]),
       boost::make_counting_iterator<const Box*>(&boxes[idx_start + polylines[i].size() - 1]));
    std::vector<const Box*> box2_ptr
      (boost::make_counting_iterator<const Box*>(&boxes[idx_start + polylines[i].size() - 1]),
       boost::make_counting_iterator<const Box*>(&boxes.back()));
    idx_start += polylines[i].size() - 1;

    CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                             box2_ptr.begin(), box2_ptr.end(),
                             intersect_polylines, cutoff);
  }
}


} } //end of namespace CGAL::Polygon_mesh_processing

/// \endcond

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H
