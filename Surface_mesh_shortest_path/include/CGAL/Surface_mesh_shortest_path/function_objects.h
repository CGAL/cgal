// Copyright (c) 2014 GeometryFactory
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
//
// Author(s)     : Stephen Kiazyk

#include <cstddef>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/result_of.h>
#include <CGAL/assertions.h>

#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H

namespace CGAL {

namespace Surface_mesh_shortest_paths_3 {

template <class Kernel>
class Compute_parametric_distance_along_segment_2
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Segment_2 Segment_2;

  typedef typename Kernel::Intersect_2 Intersect_2;
  typedef typename Kernel::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_2;
  typedef typename Kernel::Construct_vector_2 Construct_vector_2;
  typedef typename Kernel::Construct_source_2 Construct_source_2;
  typedef typename Kernel::Construct_target_2 Construct_target_2;

  typedef typename Kernel::Cartesian_const_iterator_2 Cartesian_const_iterator_2;

  typedef FT result_type;

private:
  Intersect_2 m_intersect_2;
  Construct_cartesian_const_iterator_2 m_construct_cartesian_const_iterator_2;
  Construct_vector_2 m_construct_vector_2;
  Construct_source_2 m_construct_source_2;
  Construct_target_2 m_construct_target_2;


public:
  Compute_parametric_distance_along_segment_2()
  {
  }

  Compute_parametric_distance_along_segment_2(const Kernel& kernel)
    : m_intersect_2(kernel.intersect_2_object())
    , m_construct_cartesian_const_iterator_2(kernel.construct_cartesian_const_iterator_2_object())
    , m_construct_vector_2(kernel.construct_vector_2_object())
    , m_construct_source_2(kernel.construct_source_2_object())
    , m_construct_target_2(kernel.construct_target_2_object())
  {
  }

  result_type operator () (const Segment_2& s, const Point_2& p) const
  {
    return (*this)(m_construct_source_2(s), m_construct_target_2(s), p);
  }

  result_type operator () (const Point_2& x0, const Point_2& x1, const Point_2& point) const
  {
    Vector_2 lineDiff(m_construct_vector_2(x0, x1));
    Vector_2 pointDiff(m_construct_vector_2(x0, point));

    Cartesian_const_iterator_2 lineDiffIt(m_construct_cartesian_const_iterator_2(lineDiff));
    Cartesian_const_iterator_2 pointDiffIt(m_construct_cartesian_const_iterator_2(pointDiff));

    FT lineDiff_x = *lineDiffIt;
    ++lineDiffIt;
    FT lineDiff_y = *lineDiffIt;

    FT pointDiff_x = *pointDiffIt;
    ++pointDiffIt;
    FT pointDiff_y = *pointDiffIt;

    if (CGAL::abs(lineDiff_x) > CGAL::abs(lineDiff_y))
    {
      return pointDiff_x / lineDiff_x;
    }
    else
    {
      return pointDiff_y / lineDiff_y;
    }
  }
};

template <class Kernel>
class Parametric_distance_along_segment_3
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Segment_3 Segment_3;

  typedef typename Kernel::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_3;
  typedef typename Kernel::Construct_vector_3 Construct_vector_3;
  typedef typename Kernel::Construct_source_3 Construct_source_3;
  typedef typename Kernel::Construct_target_3 Construct_target_3;

  typedef typename Kernel::Cartesian_const_iterator_3 Cartesian_const_iterator_3;

  typedef FT result_type;

private:
  Construct_cartesian_const_iterator_3 m_construct_cartesian_const_iterator_3;
  Construct_vector_3 m_construct_vector_3;
  Construct_source_3 m_construct_source_3;
  Construct_target_3 m_construct_target_3;


public:
  Parametric_distance_along_segment_3()
  {
  }

  Parametric_distance_along_segment_3(const Kernel& kernel)
    : m_construct_cartesian_const_iterator_3(kernel.construct_cartesian_const_iterator_3_object())
    , m_construct_vector_3(kernel.construct_vector_3_object())
    , m_construct_source_3(kernel.construct_source_3_object())
    , m_construct_target_3(kernel.construct_target_3_object())
  {
  }

  result_type operator () (const Segment_3& s, const Point_3& p) const
  {
    return (*this)(m_construct_source_3(s), m_construct_target_3(s), p);
  }

  result_type operator () (const Point_3& x0, const Point_3& x1, const Point_3& point) const
  {
    Vector_3 lineDiff(m_construct_vector_3(x0, x1));
    Vector_3 pointDiff(m_construct_vector_3(x0, point));

    Cartesian_const_iterator_3 lineDiffIt(m_construct_cartesian_const_iterator_3(lineDiff));
    Cartesian_const_iterator_3 pointDiffIt(m_construct_cartesian_const_iterator_3(pointDiff));

    FT lineDiff_x = *lineDiffIt;
    ++lineDiffIt;
    FT lineDiff_y = *lineDiffIt;
    ++lineDiffIt;
    FT lineDiff_z = *lineDiffIt;

    FT pointDiff_x = *pointDiffIt;
    ++pointDiffIt;
    FT pointDiff_y = *pointDiffIt;
    ++pointDiffIt;
    FT pointDiff_z = *pointDiffIt;

    if (CGAL::abs(lineDiff_x) > CGAL::abs(lineDiff_y) && CGAL::abs(lineDiff_x) > CGAL::abs(lineDiff_z))
    {
      return pointDiff_x / lineDiff_x;
    }
    else if (CGAL::abs(lineDiff_y) > CGAL::abs(lineDiff_z))
    {
      return pointDiff_y / lineDiff_y;
    }
    else
    {
      return pointDiff_z / lineDiff_z;
    }
  }
};

template<class K>
class Construct_triangle_3_to_triangle_2_projection
{
public:
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Line_3 Line_3;

  typedef typename K::Construct_line_3 Construct_line_3;
  typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;
  typedef typename K::Construct_vertex_3 Construct_vertex_3;
  typedef typename K::Construct_vector_3 Construct_vector_3;
  typedef typename K::Construct_point_2 Construct_point_2;
  typedef typename K::Construct_triangle_2 Construct_triangle_2;
  typedef typename K::Construct_projected_point_3 Construct_projected_point_3;

private:
  Parametric_distance_along_segment_3<K> m_parametric_distance_along_segment_3;
  Compute_squared_distance_3 m_compute_squared_distance_3;
  Construct_line_3 m_construct_line_3;
  Construct_projected_point_3 m_construct_projected_point_3;
  mutable Construct_vertex_3 m_construct_vertex_3;
  //~ Construct_vector_3 m_construct_vector_3;
  Construct_point_2 m_construct_point_2;
  Construct_triangle_2 m_construct_triangle_2;

public:
  Construct_triangle_3_to_triangle_2_projection()
  {
  }

  Construct_triangle_3_to_triangle_2_projection(const K& kernel)
    : m_compute_squared_distance_3(kernel.compute_squared_distance_3_object())
    , m_construct_line_3(kernel.construct_line_3_object())
    , m_construct_vertex_3(kernel.construct_vertex_3_object())
    , m_construct_projected_point_3(kernel.construct_projected_point_3_object())
    , m_construct_point_2(kernel.construct_point_2_object())
    , m_construct_triangle_2(kernel.construct_triangle_2_object())
  {
  }

  Triangle_2 operator() (const Triangle_3& t3) const
  {
    Line_3 baseSegment(m_construct_line_3(m_construct_vertex_3(t3, 0), m_construct_vertex_3(t3, 1)));

    Point_3 projectedLocation3d(m_construct_projected_point_3(baseSegment, m_construct_vertex_3(t3, 2)));
    FT scalePoint = m_parametric_distance_along_segment_3(m_construct_vertex_3(t3, 0), m_construct_vertex_3(t3, 1), projectedLocation3d);
    FT triangleHeight = CGAL::internal::select_sqrt(m_compute_squared_distance_3(projectedLocation3d, t3[2]));
    FT v01Len = CGAL::internal::select_sqrt(m_compute_squared_distance_3(t3[1], t3[0]));

    Point_2 A(m_construct_point_2(0.0, 0.0));
    Point_2 B(m_construct_point_2(v01Len, 0.0));
    Point_2 C(m_construct_point_2(v01Len * scalePoint, triangleHeight));

    return m_construct_triangle_2(A, B, C);
  }
};

template<class K>
class Robust_project_triangle_3_to_triangle_2
{
public:
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Triangle_2 Triangle_2;
  typedef Exact_predicates_exact_constructions_kernel_with_sqrt EKSQRT;
  typedef Construct_triangle_3_to_triangle_2_projection<EKSQRT> Exact_project_triangle_3_to_triangle_2;
  typedef Cartesian_converter<K, EKSQRT>  To_exact;
  typedef Cartesian_converter<EKSQRT, K>  Back_from_exact;

public:
  Robust_project_triangle_3_to_triangle_2()
  {
  }

  Robust_project_triangle_3_to_triangle_2(const K& /* kernel */)
  {
  }

  Triangle_2 operator() (const Triangle_3& t3) const
  {
    Exact_project_triangle_3_to_triangle_2 ept3t2;
    To_exact to_exact;
    Back_from_exact back_from_exact;

    return back_from_exact(ept3t2(to_exact(t3)));
  }
};

template<class K>
class Construct_triangle_3_along_segment_2_flattening
{
public:
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Line_3 Line_3;
  typedef typename K::FT FT;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_2 Segment_2;

  typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;
  typedef typename K::Construct_projected_point_3 Construct_projected_point_3;
  typedef typename K::Construct_perpendicular_vector_2 Construct_perpendicular_vector_2;
  typedef typename K::Construct_sum_of_vectors_2 Construct_sum_of_vectors_2;
  typedef typename K::Construct_scaled_vector_2 Construct_scaled_vector_2;
  typedef typename K::Construct_translated_point_2 Construct_translated_point_2;
  typedef typename K::Construct_vector_2 Construct_vector_2;
  typedef typename K::Compute_squared_length_2 Compute_squared_length_2;
  typedef typename K::Construct_vertex_3 Construct_vertex_3;
  typedef typename K::Construct_triangle_2 Construct_triangle_2;
  typedef typename K::Construct_line_3 Construct_line_3;
  typedef typename K::Construct_segment_3 Construct_segment_3;
  typedef typename K::Construct_source_2 Construct_source_2;
  typedef typename K::Construct_target_2 Construct_target_2;

  typedef Triangle_2 result_type;

private:

  Parametric_distance_along_segment_3<K> m_parametric_distance_along_segment_3;
  Compute_squared_distance_3 m_compute_squared_distance_3;
  Construct_projected_point_3 m_construct_projected_point_3;
  Construct_perpendicular_vector_2 m_construct_perpendicular_vector_2;
  Construct_sum_of_vectors_2 m_construct_sum_of_vectors_2;
  Construct_scaled_vector_2 m_construct_scaled_vector_2;
  Construct_translated_point_2 m_construct_translated_point_2;
  Construct_vector_2 m_construct_vector_2;
  Compute_squared_length_2 m_compute_squared_length_2;
  Construct_line_3 m_construct_line_3;
  Construct_segment_3 m_construct_segment_3;
  Construct_source_2 m_construct_source_2;
  Construct_target_2 m_construct_target_2;
  mutable Construct_vertex_3 m_construct_vertex_3;
  Construct_triangle_2 m_construct_triangle_2;

public:
  Construct_triangle_3_along_segment_2_flattening()
  {
  }

  Construct_triangle_3_along_segment_2_flattening(const K& kernel)
    : m_compute_squared_distance_3(kernel.compute_squared_distance_3_object())
    , m_construct_projected_point_3(kernel.construct_projected_point_3_object())
    , m_construct_perpendicular_vector_2(kernel.construct_perpendicular_vector_2_object())
    , m_construct_sum_of_vectors_2(kernel.construct_sum_of_vectors_2_object())
    , m_construct_scaled_vector_2(kernel.construct_scaled_vector_2_object())
    , m_construct_translated_point_2(kernel.construct_translated_point_2_object())
    , m_construct_vector_2(kernel.construct_vector_2_object())
    , m_compute_squared_length_2(kernel.compute_squared_length_2_object())
    , m_construct_line_3(kernel.construct_line_3_object())
    , m_construct_segment_3(kernel.construct_segment_3_object())
    , m_construct_source_2(kernel.construct_source_2_object())
    , m_construct_target_2(kernel.construct_target_2_object())
    , m_construct_vertex_3(kernel.construct_vertex_3_object())
    , m_construct_triangle_2(kernel.construct_triangle_2_object())
  {
  }

  result_type operator() (const Triangle_3& t3, std::size_t edgeIndex, const Segment_2& segment) const
  {
    Point_3 projectedLocation3d(m_construct_projected_point_3(m_construct_line_3(m_construct_vertex_3(t3, edgeIndex), m_construct_vertex_3(t3, edgeIndex + 1)), m_construct_vertex_3(t3, edgeIndex + 2)));
    FT scalePoint = m_parametric_distance_along_segment_3(m_construct_segment_3(m_construct_vertex_3(t3, edgeIndex), m_construct_vertex_3(t3, edgeIndex + 1)), projectedLocation3d);
    FT triangleHeight = CGAL::internal::select_sqrt(m_compute_squared_distance_3(projectedLocation3d, m_construct_vertex_3(t3, edgeIndex + 2)));

    Vector_2 edgeVector(m_construct_vector_2(segment));

    Vector_2 perpendicularEdgeVector(m_construct_perpendicular_vector_2(edgeVector, CGAL::COUNTERCLOCKWISE));
    perpendicularEdgeVector = m_construct_scaled_vector_2(perpendicularEdgeVector, FT(1.0) / CGAL::internal::select_sqrt(m_compute_squared_length_2(perpendicularEdgeVector)));

    Point_2 points[3];
    points[edgeIndex] = m_construct_source_2(segment);
    points[(edgeIndex + 1) % 3] = m_construct_target_2(segment);
    points[(edgeIndex + 2) % 3] = m_construct_translated_point_2(m_construct_source_2(segment), m_construct_sum_of_vectors_2(m_construct_scaled_vector_2(edgeVector, scalePoint), m_construct_scaled_vector_2(perpendicularEdgeVector, triangleHeight)));

    return m_construct_triangle_2(points[0], points[1], points[2]);
  }
};

template<class K>
class Robust_flatten_triangle_3_along_segment_2
{
public:
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Triangle_2 Triangle_2;

  typedef Exact_predicates_exact_constructions_kernel_with_sqrt EKSQRT;
  typedef Construct_triangle_3_along_segment_2_flattening<EKSQRT> Exact_flatten_triangle_3_along_segment_2;
  typedef Cartesian_converter<K, EKSQRT>  To_exact;
  typedef Cartesian_converter<EKSQRT, K>  Back_from_exact;

public:
  Robust_flatten_triangle_3_along_segment_2()
  {
  }

  Robust_flatten_triangle_3_along_segment_2(const K& /* kernel */)
  {
  }

  Triangle_2 operator() (const Triangle_3& t3, std::size_t edgeIndex, const Segment_2& segment) const
  {
    Exact_flatten_triangle_3_along_segment_2 eft3as2;
    To_exact to_exact;
    Back_from_exact back_from_exact;

    return back_from_exact(eft3as2(to_exact(t3), edgeIndex, to_exact(segment)));
  }
};

template <class K>
class Compare_relative_intersection_along_segment_2
{
public:
  typedef typename K::FT FT;
  typedef typename K::Ray_2 Ray_2;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Line_2 Line_2;

  typedef typename K::Intersect_2 Intersect_2;
  typedef typename K::Compare_distance_2 Compare_distance_2;
  typedef typename K::Construct_line_2 Construct_line_2;
  typedef typename K::Construct_source_2 Construct_source_2;
  typedef typename K::Construct_target_2 Construct_target_2;

  typedef CGAL::Comparison_result result_type;

private:
  Compute_parametric_distance_along_segment_2<K> m_parametric_distance_along_segment_2;
  Intersect_2 m_intersect_2;
  Compare_distance_2 m_compare_distance_2;
  Construct_line_2 m_construct_line_2;
  Construct_source_2 m_construct_source_2;
  Construct_target_2 m_construct_target_2;

public:
  Compare_relative_intersection_along_segment_2()
  {
  }

  Compare_relative_intersection_along_segment_2(const K& kernel)
    : m_intersect_2(kernel.intersect_2_object())
    , m_compare_distance_2(kernel.compare_distance_2_object())
    , m_construct_line_2(kernel.construct_line_2_object())
    , m_construct_source_2(kernel.construct_source_2_object())
    , m_construct_target_2(kernel.construct_target_2_object())
  {
  }

  result_type operator () (const Segment_2& s1, const Line_2& l1, const Segment_2& s2, const Line_2& l2) const
  {
    typedef typename CGAL::cpp11::result_of<Intersect_2(Line_2, Line_2)>::type LineLineIntersectResult;

    Line_2 s1Line(m_construct_line_2(s1));

    Line_2 s2Line(m_construct_line_2(s2));

    LineLineIntersectResult intersectResult1(m_intersect_2(s1Line, l1));

    CGAL_assertion(bool(intersectResult1));

    Point_2 p1;

    FT t1;

    if (intersectResult1)
    {
      Point_2* result = boost::get<Point_2, Point_2, Line_2>(&*intersectResult1);

      CGAL_assertion(result && "Intersection should have been a point");

      if (result)
      {
        t1 = m_parametric_distance_along_segment_2(s1, *result);
        p1 = *result;
        CGAL_assertion(t1 >= FT(-0.00001) && t1 <= FT(1.00001));
      }
    }

    LineLineIntersectResult intersectResult2 = m_intersect_2(s2Line, l2);

    CGAL_assertion(bool(intersectResult2));

    FT t2;
    Point_2 p2;

    if (intersectResult2)
    {
      Point_2* result = boost::get<Point_2, Point_2, Line_2>(&*intersectResult2);

      CGAL_assertion(result && "Intersection should have been a point");

      if (result)
      {
        t2 = m_parametric_distance_along_segment_2(s2, *result);
        p2 = *result;
        CGAL_assertion(t2 >= FT(-0.00001) && t2 <= FT(1.00001));
      }
    }

    result_type predicateResult = m_compare_distance_2(s1.source(), p1, s2.source(), p2);

    return predicateResult;
  }
};

template <class Kernel, class FaceListGraph>
class Is_saddle_vertex
{
public:
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_2 Point_2;

  typedef typename boost::graph_traits<FaceListGraph> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  typedef typename CGAL::Surface_mesh_shortest_paths_3::Construct_triangle_3_to_triangle_2_projection<Kernel> Construct_triangle_3_to_triangle_2_projection;
  typedef typename CGAL::Surface_mesh_shortest_paths_3::Construct_triangle_3_along_segment_2_flattening<Kernel> Construct_triangle_3_along_segment_2_flattening;
  typedef typename Kernel::Orientation_2 Orientation_2;
  typedef typename Kernel::Construct_triangle_3 Construct_triangle_3;
  typedef typename Kernel::Construct_vertex_2 Construct_vertex_2;
  typedef typename Kernel::Construct_segment_2 Construct_segment_2;
  typedef typename Kernel::Construct_source_2 Construct_source_2;
  typedef typename Kernel::Construct_target_2 Construct_target_2;

  typedef typename Kernel::Boolean result_type;

private:
  Construct_triangle_3_to_triangle_2_projection m_project_triangle_3_to_triangle_2;
  Construct_triangle_3_along_segment_2_flattening m_flatten_triangle_3_along_segment_2;
  Construct_triangle_3 m_construct_triangle_3;
  mutable Construct_vertex_2 m_construct_vertex_2;
  Construct_segment_2 m_construct_segment_2;
  Construct_source_2 m_construct_source_2;
  Construct_target_2 m_construct_target_2;
  Orientation_2 m_orientation_2;

public:

  Is_saddle_vertex()
  {
  }

  Is_saddle_vertex(const Kernel& kernel, const Construct_triangle_3_to_triangle_2_projection& pt3tt2, const Construct_triangle_3_along_segment_2_flattening& ft3as2)
    : m_orientation_2(kernel.orientation_2_object())
    , m_construct_triangle_3(kernel.construct_triangle_3_object())
    , m_construct_vertex_2(kernel.construct_vertex_2_object())
    , m_construct_segment_2(kernel.construct_segment_2_object())
    , m_construct_source_2(kernel.construct_source_2_object())
    , m_construct_target_2(kernel.construct_target_2_object())
    , m_project_triangle_3_to_triangle_2(pt3tt2)
    , m_flatten_triangle_3_along_segment_2(ft3as2)
  {
  }

  result_type operator() (vertex_descriptor v, FaceListGraph& g) const
  {
    return (*this)(v, g, get(boost::vertex_point, g));
  }

  template<class VertexPointMap>
  result_type operator() (vertex_descriptor v, const FaceListGraph& g, VertexPointMap const& pointMap) const
  {
    halfedge_descriptor startEdge = halfedge(v, g);

    Point_3 rootPoint(get(pointMap, v));
    Point_3 prevPoint(get(pointMap, source(startEdge, g)));

    halfedge_descriptor currentEdge = next(startEdge, g);

    Point_3 nextPoint(get(pointMap, target(currentEdge, g)));

    Triangle_3 baseFace3(rootPoint, nextPoint, prevPoint);

    currentEdge = opposite(currentEdge, g);

    Triangle_2 baseFace2(m_project_triangle_3_to_triangle_2(baseFace3));

    Point_2 A(m_construct_vertex_2(baseFace2, 0));
    Point_2 B(m_construct_vertex_2(baseFace2, 1));
    Point_2 C(m_construct_vertex_2(baseFace2, 2));

    Segment_2 baseSegment(m_construct_segment_2(A, C));

    Segment_2 nextSegment(m_construct_segment_2(B, A));

    CGAL::Orientation baseOrientation = m_orientation_2(m_construct_vertex_2(baseFace2, 0), m_construct_vertex_2(baseFace2, 2), m_construct_vertex_2(baseFace2, 1));

    CGAL_assertion(baseOrientation != CGAL::COLLINEAR);

    do
    {
      prevPoint = nextPoint;
      currentEdge = next(currentEdge, g);
      nextPoint = get(pointMap, target(currentEdge, g));
      currentEdge = opposite(currentEdge, g);

      Triangle_3 currentFace3(m_construct_triangle_3(rootPoint, nextPoint, prevPoint));
      Triangle_2 currentFace2(m_flatten_triangle_3_along_segment_2(currentFace3, 2, nextSegment));

      if (m_orientation_2(m_construct_source_2(baseSegment), m_construct_target_2(baseSegment), m_construct_vertex_2(currentFace2, 2)) != baseOrientation && m_orientation_2(m_construct_source_2(baseSegment), m_construct_target_2(baseSegment), m_construct_vertex_2(currentFace2, 1)) == baseOrientation)
      {
        return true;
      }

      nextSegment = m_construct_segment_2(currentFace2[1], currentFace2[0]);
    }
    while (currentEdge != startEdge);

    return false;
  }
};

} // namespace Surface_mesh_shortest_paths_3

} // namespace CGAL

#endif /* CGAL_SURFACE_MESH_SHORTEST_PATHS_3_FUNCTION_OBJECTS_H */
