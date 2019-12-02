// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H

#include <CGAL/license/Surface_mesh_shortest_path.h>

#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#ifdef CGAL_SMSP_USE_ROBUST_TRAITS_CODE
#if defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#endif
#endif

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>
#include <CGAL/result_of.h>
#include <CGAL/Cartesian_converter.h>

#include <cmath>
#include <cstddef>
#include <limits>

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
    FT triangleHeight = CGAL::approximate_sqrt(m_compute_squared_distance_3(projectedLocation3d, t3[2]));
    FT v01Len = CGAL::approximate_sqrt(m_compute_squared_distance_3(t3[1], t3[0]));

    Point_2 A(m_construct_point_2(0.0, 0.0));
    Point_2 B(m_construct_point_2(v01Len, 0.0));
    Point_2 C(m_construct_point_2(v01Len * scalePoint, triangleHeight));

    return m_construct_triangle_2(A, B, C);
  }
};

#ifdef CGAL_SMSP_USE_ROBUST_TRAITS_CODE
#if defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)
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
#endif
#endif

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

  result_type operator() (const Triangle_3& t3, int edgeIndex, const Segment_2& segment) const
  {
    Point_3 projectedLocation3d(m_construct_projected_point_3(m_construct_line_3(m_construct_vertex_3(t3, edgeIndex), m_construct_vertex_3(t3, edgeIndex + 1)), m_construct_vertex_3(t3, edgeIndex + 2)));
    FT scalePoint = m_parametric_distance_along_segment_3(m_construct_segment_3(m_construct_vertex_3(t3, edgeIndex), m_construct_vertex_3(t3, edgeIndex + 1)), projectedLocation3d);
    FT triangleHeight = CGAL::approximate_sqrt(m_compute_squared_distance_3(projectedLocation3d, m_construct_vertex_3(t3, edgeIndex + 2)));

    Vector_2 edgeVector(m_construct_vector_2(segment));

    Vector_2 perpendicularEdgeVector(m_construct_perpendicular_vector_2(edgeVector, CGAL::COUNTERCLOCKWISE));
    perpendicularEdgeVector = m_construct_scaled_vector_2(perpendicularEdgeVector, FT(1) / CGAL::approximate_sqrt(m_compute_squared_length_2(perpendicularEdgeVector)));

    Point_2 points[3];
    points[edgeIndex] = m_construct_source_2(segment);
    points[(edgeIndex + 1) % 3] = m_construct_target_2(segment);
    points[(edgeIndex + 2) % 3] = m_construct_translated_point_2(m_construct_source_2(segment), m_construct_sum_of_vectors_2(m_construct_scaled_vector_2(edgeVector, scalePoint), m_construct_scaled_vector_2(perpendicularEdgeVector, triangleHeight)));

    return m_construct_triangle_2(points[0], points[1], points[2]);
  }
};

#ifdef CGAL_SMSP_USE_ROBUST_TRAITS_CODE
#if defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)
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
#endif
#endif

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
  typedef typename K::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef typename K::Construct_line_2 Construct_line_2;
  typedef typename K::Construct_source_2 Construct_source_2;
  typedef typename K::Construct_target_2 Construct_target_2;


  typedef CGAL::Comparison_result result_type;

private:
  Compute_parametric_distance_along_segment_2<K> m_parametric_distance_along_segment_2;
  Intersect_2 m_intersect_2;
  Compare_distance_2 m_compare_distance_2;
  Compute_squared_distance_2 m_compute_squared_distance_2;
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
    , m_compute_squared_distance_2(kernel.compute_squared_distance_2_object())
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
    if (!intersectResult1) return CGAL::SMALLER;

    const Point_2* p1_ptr = boost::get<Point_2>(&*intersectResult1);

    CGAL_assertion(p1_ptr && "Intersection should have been a point");
    if (!p1_ptr) return CGAL::SMALLER;

    CGAL_assertion_code(FT t1 = m_parametric_distance_along_segment_2(s1, *p1_ptr);)
    CGAL_assertion(t1 >= FT(-1)/FT(100000) && t1 <= FT(1)+FT(1)/FT(100000));

    LineLineIntersectResult intersectResult2 = m_intersect_2(s2Line, l2);
    CGAL_assertion(bool(intersectResult2));
    if (!intersectResult2) return CGAL::SMALLER;

    const Point_2* p2_ptr = boost::get<Point_2>(&*intersectResult2);

    CGAL_assertion(p2_ptr && "Intersection should have been a point");
    if (!p2_ptr) return CGAL::SMALLER;

    CGAL_assertion_code(FT t2 = m_parametric_distance_along_segment_2(s2, *p2_ptr);)
    CGAL_assertion(t2 >= FT(-1)/FT(100000) && t2 <= FT(1)+FT(1)/FT(100000));

// #define CGAL_SMSP_DONT_USE_RELAXED_PRUNING
#ifndef CGAL_SMSP_DONT_USE_RELAXED_PRUNING
    const FT sqd_1 = m_compute_squared_distance_2(s1.source(), *p1_ptr);
    const FT sqd_2 = m_compute_squared_distance_2(s2.source(), *p2_ptr);

    // In the case of multiple rays reaching the same target, we want to know their respective position
    // so that pruning of branches can be done according to the "one angle one split" idiom.
    // However, the orientation predicate is evaluated in the unfolded 2D plane, which is obtained
    // via square roots; inconsisnties will exist. We don't want to prune in case it might be wrong,
    // so we add a little bit of tolerance on the evaluation of the predicate. If it's almost collinear,
    // return 'collinear' (EQUAL).
    const FT eps = (FT(100) * std::numeric_limits<FT>::epsilon());
    if(CGAL::abs(sqd_1 - sqd_2) < eps)
      return CGAL::EQUAL;

    return CGAL::compare(sqd_1, sqd_2);
#else
    return m_compare_distance_2(s1.source(), *p1_ptr, s2.source(), *p2_ptr);
#endif
  }
};

template <class Kernel, class FaceListGraph>
class Is_saddle_vertex
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Triangle_2 Triangle_2;

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

  Is_saddle_vertex(const Kernel& kernel,
                   const Construct_triangle_3_to_triangle_2_projection& pt3tt2,
                   const Construct_triangle_3_along_segment_2_flattening& ft3as2)
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

 FT angle(const Vector_3& u, const Vector_3& v) const
 {
   typename Kernel::Compute_scalar_product_3 scalar_product;

   double product = CGAL::sqrt(to_double(scalar_product(u,u)) * to_double(scalar_product(v,v)));

   if(product == 0)
     return 0;

   // cosine
   double dot = to_double(scalar_product(u,v));
   double cosine = dot / product;

   if(cosine > 1.)
     cosine = 1.;

   if(cosine < -1.)
     cosine = -1.;

   return std::acos(cosine) * 180./CGAL_PI;
 }


 FT angle(const Point_3& p, const Point_3& q, const Point_3& r) const
 {
   typename Kernel::Construct_vector_3 cv;

   Vector_3 u = cv(q, p);
   Vector_3 v = cv(q, r);

   return angle(u, v);
 }

  template<class VertexPointMap>
  FT vertex_angle(const vertex_descriptor v,
                  const FaceListGraph& g,
                  const VertexPointMap& pointMap) const
  {
    FT angle_sum = 0;

    for(halfedge_descriptor h : halfedges_around_target(v, g))
    {
      if(is_border(h, g))
        continue;

      angle_sum += angle(get(pointMap, source(h, g)),
                         get(pointMap, target(h, g)),
                         get(pointMap, target(next(h, g), g)));
    }

    angle_sum *= CGAL_PI / FT(180);

    return angle_sum;
  }

  template<class VertexPointMap>
  result_type operator() (vertex_descriptor v, const FaceListGraph& g, VertexPointMap const& pointMap) const
  {
#ifndef CGAL_SMSP_DONT_USE_RELAXED_PRUNING
    const FT ang_sum = vertex_angle(v, g, pointMap);
    const FT bound = (FT(1) - FT(100) * std::numeric_limits<FT>::epsilon()) * 2 * CGAL_PI;
    return (ang_sum >= bound);
#else
    halfedge_descriptor startEdge = halfedge(v, g);
    while (face(startEdge, g) == Graph_traits::null_face())
      startEdge=opposite(next(startEdge, g), g);

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
#endif
  }
};

} // namespace Surface_mesh_shortest_paths_3

} // namespace CGAL

#endif /* CGAL_SURFACE_MESH_SHORTEST_PATHS_3_FUNCTION_OBJECTS_H */
