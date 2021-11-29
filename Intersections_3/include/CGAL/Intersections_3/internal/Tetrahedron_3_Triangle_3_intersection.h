// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H

#include <CGAL/Intersections_3/internal/Plane_3_Tetrahedron_3_intersection.h>
#include <CGAL/Intersections_3/internal/Triangle_3_Triangle_3_intersection.h>

#include <CGAL/kernel_basic.h>

#include <algorithm>
#include <iterator>
#include <list>
#include <utility>
#include <vector>

namespace CGAL {
namespace Intersections {
namespace internal {

template<typename Segment>
void filter_segments(std::list<Segment>& segments)
{
  auto are_equal = [](const Segment& l, const Segment& r) -> bool
                   {
                     return (l == r || l == r.opposite());
                   };

  auto it = std::unique(segments.begin(), segments.end(), are_equal);
  segments.erase(it, segments.end());
}

// plane going through the ref segment's source, a point above (given by the normal of the input
// triangle) and two points (ref_other / query) that need to be ordered
template <class K>
bool first_comes_first_pt(const typename K::Point_3& ref_source,
                          const typename K::Point_3& ref_z,
                          const typename K::Point_3& ref_other,
                          const typename K::Point_3& query,
                          const K& k)
{
  typename K::Orientation_3 orientation = k.orientation_3_object();

  // points have filtered to remove segments' extremities
  CGAL_precondition(ref_other != query);

  const Orientation o = orientation(ref_source, ref_z, ref_other, query);
  CGAL_assertion(o != COPLANAR);

  // ref_other comes first <==> query is on the positive side of the plane
  return (o == POSITIVE);
}

template <class K, class SegPtVariant>
bool first_comes_first(const typename K::Point_3& ref_source,
                       const typename K::Point_3& ref_z,
                       const typename K::Point_3& ref_other,
                       const SegPtVariant& seg_or_pt,
                       const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;

  typedef typename std::list<Segment_3>::iterator SCI;
  typedef typename std::vector<Point_3>::iterator PCI;

  if(seg_or_pt.which() == 0)
  {
    const Segment_3& s = *(boost::get<SCI>(seg_or_pt));
    return first_comes_first_pt(ref_source, ref_z, ref_other, s.source(), k);
  }
  else
  {
    CGAL_assertion(seg_or_pt.which() == 1);

    const Point_3& p = *(boost::get<PCI>(seg_or_pt));
    return first_comes_first_pt(ref_source, ref_z, ref_other, p, k);
  }
}

template <class K, class SegmentContainer, class PointContainer>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
build_intersection(const typename K::Tetrahedron_3& /*input_tetrahedron*/,
                   const typename K::Triangle_3& input_triangle,
                   PointContainer& points,
                   SegmentContainer& segments,
                   const K& k)
{
  typedef typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type result_type;

  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef std::vector<Point_3> Poly;

  typedef typename SegmentContainer::iterator SCI;
  typedef typename PointContainer::iterator PCI;

  // @todo? Could do the 1 segment case with this code too...
  CGAL_precondition(segments.size() >= 2 && segments.size() <= 4);
  CGAL_precondition(points.size() <= 2);

  // Constructions @fixme avoidable?
  const Vector_3 input_triangle_normal = input_triangle.supporting_plane().orthogonal_vector();

  // remove points that are just segments extremities
  auto is_extremity = [&segments](const Point_3& p) -> bool
                      {
                        for(const Segment_3& s : segments)
                          if(p == s.source() || p == s.target())
                            return true;
                        return false;
                      };
  points.erase(std::remove_if(points.begin(), points.end(), is_extremity),
               points.end());

  // Take the first segment as reference, and order the rest to form a convex polygon
  //
  // All segments and points involved in the intersection are on the input triangle
  // and thus everything is coplanar, at least theoretically (the kernel might not provide
  // exact constructions...)
  //
  // Given an arbitrary segment, the code below sorts the other segments and the points
  // in a ccw order. Using a vector because the number of segments and points is bounded
  // (max 4 segments and max 2 points) so even if the linear insertion is a little ugly,
  // it is not expensive anyway.
  //
  // Example:
  /*
               x p0

           /
          /
      s1 /          \ s2
        /            \
        --------------
              s0
  */
  //
  // s0 is chosen as the reference segment
  // output will be 's0 s2 p0 s1'

  Segment_3& ref_s = segments.front();
  Point_3 ref_z = ref_s.source() + input_triangle_normal;

  // The reference segment should be such that all other intersection parts are
  // on the positive side of the plane described by the normal of the triangle and ref_s
  bool swapped = false;
  for(SCI slit = std::next(segments.begin()); slit != segments.end(); ++slit)
  {
    const Segment_3& other = *slit;

    if(k.orientation_3_object()(ref_s.source(), ref_z, ref_s.target(), other.source()) == CGAL::NEGATIVE ||
       k.orientation_3_object()(ref_s.source(), ref_z, ref_s.target(), other.target()) == CGAL::NEGATIVE)
    {
      ref_s = ref_s.opposite();
      ref_z = ref_s.source() + input_triangle_normal;
      swapped = true;
      break;
    }
  }

  if(!swapped)
  {
    for(PCI plit = points.begin(); plit != points.end(); ++plit)
    {
      const Point_3& other = *plit;
      if(k.orientation_3_object()(ref_s.source(), ref_z, ref_s.target(), other) == CGAL::NEGATIVE)
      {
        swapped = true;
        ref_s = ref_s.opposite();
        ref_z = ref_s.source() + input_triangle_normal;
        break;
      }
    }
  }

  const Point_3& ref_sp = ref_s.source();
  const Point_3& ref_tp = ref_s.target();

  // Now, order the other parts of the intersection
  std::list<boost::variant<SCI, PCI> > res_elements; // iterators to the points/segments
  res_elements.emplace_back(segments.begin());

  for(SCI slit = std::next(segments.begin()); slit != segments.end(); ++slit)
  {
    // first, check if the segment is well oriented, meaning its source comes before its target (ccwly)
    Segment_3& curr_s = *slit;

    if(curr_s.source() == ref_sp || curr_s.target() == ref_tp) // consecutive segments must have consistent orientation
    {
      curr_s = curr_s.opposite();
    }
    else if(curr_s.source() == ref_tp || curr_s.target() == ref_sp)
    {
      // nothing to do here as we know that sp&tp are on the positive side of (normal, ref_s)
    }
    else if(first_comes_first_pt(ref_sp, ref_z, curr_s.target(), curr_s.source(), k))
    {
      curr_s = curr_s.opposite();
    }

    // Find where the current segment fit in the final polygon intersection
    for(auto rit = std::next(res_elements.begin()); ; ++rit)
    {
      // always pick the current segment's source to ensure ref_source != ref_other
      if(rit == res_elements.end() || first_comes_first(ref_sp, ref_z, curr_s.source(), *rit, k))
      {
        res_elements.insert(rit, slit);
        break;
      }
    }
  }

  for(PCI plit = points.begin(); plit != points.end(); ++plit)
  {
    const Point_3& curr_p = *plit;

    // Find where the current point fits in the boundary of the polygon intersection
    for(auto rit = std::next(res_elements.begin()); ; ++rit)
    {
      if(rit == res_elements.end() || first_comes_first(ref_sp, ref_z, curr_p, *rit, k))
      {
        res_elements.insert(rit, plit);
        break;
      }
    }
  }

  CGAL_postcondition(res_elements.size() == points.size() + segments.size());

  // Concatenate points to create the polygonal output
  Poly res;
  for(const boost::variant<SCI, PCI>& e : res_elements)
  {
    if(const SCI* sci = boost::get<SCI>(&e))
    {
      const Segment_3& s = **sci;

      if(res.empty() || s.source() != res.back()) // common extremity for consecutive segments
        res.push_back(s.source());
      if(res.empty() || s.target() != res.front())
        res.push_back(s.target());
    }
    else if(const PCI* pci = boost::get<PCI>(&e))
    {
      res.push_back(**pci);
    }
    else
    {
      CGAL_assertion(false);
    }
  }

  CGAL_assertion(std::set<Point_3>(res.begin(), res.end()).size() == res.size());
  CGAL_assertion(res.size() >= 3);

  if(res.size() == 3)
  {
    Triangle_3 tr { res[0], res[1], res[2] };
    return result_type(std::forward<Triangle_3>(tr));
  }
  else
  {
    return result_type(std::forward<Poly>(res));
  }
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(const typename K::Tetrahedron_3& tet,
             const typename K::Triangle_3& tr,
             const K& k)
{
  typedef typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type result_type;
  typedef typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3>::result_type Inter_type;

  CGAL_precondition(!tet.is_degenerate());
  CGAL_precondition(!tr.is_degenerate());

  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef std::vector<Point_3> Poly;

  typename K::Bounded_side_3 bounded_side = k.bounded_side_3_object();
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Construct_triangle_3 triangle = k.construct_triangle_3_object();

  std::vector<Bounded_side> vertex_sides(3);

  std::vector<Point_3> points;
  int inside_points = 0;
  int strictly_inside_points = 0;

  for(int i=0; i<3; ++i)
  {
    vertex_sides[i] = bounded_side(tet, vertex(tr, i));

    if(vertex_sides[i] != ON_UNBOUNDED_SIDE)
      ++inside_points;

    if(vertex_sides[i] == ON_BOUNDED_SIDE)
    {
      ++strictly_inside_points;
      points.push_back(vertex(tr, i));
    }
  }

  switch(inside_points)
  {
    case 0:
    {
      Inter_type intersections[4];
      std::list<Segment_3> segments;
      std::vector<std::size_t> seg_ids;
      for(std::size_t i = 0; i < 4; ++i)
      {
        const Triangle_3 face = triangle(vertex(tet, (i+1)%4),
                                         vertex(tet, (i+2)%4),
                                         vertex(tet, (i+3)%4));
        intersections[i] = intersection(tr, face, k);
        if(intersections[i])
        {
          // a face is inside the input tr
          if(const Triangle_3* t = boost::get<Triangle_3>(&*intersections[i]))
          {
            Triangle_3 res = *t;
            return result_type(std::forward<Triangle_3>(res));
          }
          else if(const Segment_3* s = boost::get<Segment_3>(&*intersections[i]))
          {
            // get segs and pts to construct poly
            segments.push_back(*s);
            seg_ids.push_back(i);
          }
          else if(const Point_3* p = boost::get<Point_3>(&*intersections[i]))
          {
            points.push_back(*p);
          }
          else if(const Poly* p = boost::get<Poly>(&*intersections[i]))
          {
            // the input triangle is in the supporting plane of a tet face, return the poly.
            Poly res = *p;
            return result_type(std::forward<Poly>(res));
          }
        }
      }

      if(segments.size() > 1)
        filter_segments(segments);

      // no segments and no inside points, there can still be an intersecting (tet vertex on
      // an edge|face of the triangle)
      if(segments.empty())
      {
        if(points.empty())
          return result_type();

        return result_type(std::forward<Point_3>(points.front()));
      }
      else if(segments.size() == 1)
      {
        // adjacency to an edge, return resulting segment.
        return result_type(segments.front());
      }
      else if(segments.size() > 1)
      {
        return build_intersection(tet, tr, points, segments, k);
      }
    }
      break;
    case 1:
    case 2: // 1 or 2 inside points
    {
      Inter_type intersections[4];
      std::list<Segment_3> segments;
      for(std::size_t i = 0; i < 4; ++i)
      {
        const Triangle_3 face = triangle(vertex(tet, (i+1)%4),
                                         vertex(tet, (i+2)%4),
                                         vertex(tet, (i+3)%4));
        intersections[i] = intersection(tr, face, k);
        if(intersections[i])
        {
          if(const Triangle_3* t = boost::get<Triangle_3>(&*intersections[i]))
          {
            Triangle_3 res = *t;
            return result_type(std::forward<Triangle_3>(res));
          }
          else if(const Segment_3* s = boost::get<Segment_3>(&*intersections[i]))
          {
            segments.push_back(*s);
          }
          else if(const Point_3* p = boost::get<Point_3>(&*intersections[i]))
          {
            points.push_back(*p);
          }
          else if(const Poly* p = boost::get<Poly>(&*intersections[i]))
          {
            // the input is in a supporting plane of a face
            Poly res = *p;
            return result_type(std::forward<Poly>(res));
          }
        }
      }

      if(segments.size() > 1)
        filter_segments(segments);

      switch(segments.size())
      {
        case 0:
        {
          // there can only be one point of contact, otherwise by convexity
          // there would be a full segment on a face (interior segment isn't possible either
          // because there are at most 2 inside points and an interior segment would also
          // yield at least a segment on the boundary)
          return result_type(std::forward<Point_3>(points.front()));
        }
        case 1: // 1 segment
        {
          const Segment_3& s = segments.front();

          if(strictly_inside_points == 1)
          {
            // Doesn't matter whether there is another (non-strictly) inside point: if there is,
            // it is an extremity of the segment

            const int str_inside_pt_pos =
                int(std::find(vertex_sides.begin(), vertex_sides.end(), ON_BOUNDED_SIDE) - vertex_sides.begin());
            CGAL_assertion(str_inside_pt_pos >= 0 && str_inside_pt_pos < 3);

            Triangle_3 res_tr = triangle(vertex(tr, str_inside_pt_pos), s.source(), s.target());
            return result_type(std::forward<Triangle_3>(res_tr));
          }
          else if(strictly_inside_points == 2)
          {
            CGAL_assertion(inside_points == 2); // can't be 3 since we're in the 1&2 switch

            Poly res(4);

            // Grab the 2 strictly inside points
            int id = 0;
            for(int i=0; i<3; ++i)
              if(vertex_sides[i] == ON_BOUNDED_SIDE)
                res[id++] = vertex(tr, i);

            CGAL_assertion(id == 2);

            if((res[0] - res[1]) * (s.source() - s.target()) > 0)
            {
              res[2] = s.target();
              res[3] = s.source();
            }
            else
            {
              res[3] = s.target();
              res[2] = s.source();
            }

            return result_type(std::forward<Poly>(res));
          }
          else if(inside_points == 1) // 1 point on the boundary
          {
            CGAL_assertion(strictly_inside_points == 0);

            // Grab the inside point
            const int boundary_pt_pos =
                int(std::find(vertex_sides.begin(), vertex_sides.end(), ON_BOUNDARY) - vertex_sides.begin());
            CGAL_assertion(boundary_pt_pos >= 0 && boundary_pt_pos < 3);

            const Point_3& boundary_pt = vertex(tr, boundary_pt_pos);
            if(boundary_pt == s.source() || boundary_pt == s.target())
            {
              return result_type(s);
            }
            else
            {
              Triangle_3 res_tr = triangle(boundary_pt, s.source(), s.target());
              return result_type(std::forward<Triangle_3>(res_tr));
            }
          }
          else // 2 points on the boundary
          {
            CGAL_assertion(inside_points == 2 && strictly_inside_points == 0);

            // 2 boundary points and 1 segment, have to distinguish between cases
            // depending on if the extremities of the segment are triangle extremities

            std::array<int, 2> boundary_pts;
            std::array<bool, 2> is_boundary_point_an_extremity;

            // Grab the inside points
            std::size_t id = 0;
            for(int i=0; i<3; ++i)
            {
              if(vertex_sides[i] == ON_BOUNDARY)
              {
                boundary_pts[id] = i;

                if(vertex(tr, i) == s.source())
                  is_boundary_point_an_extremity[id] = true;
                else if(vertex(tr, i) == s.target())
                  is_boundary_point_an_extremity[id] = true;
                else
                  is_boundary_point_an_extremity[id] = false;

                ++id;
              }
            }

            CGAL_assertion(id == 2);

            if(is_boundary_point_an_extremity[0])
            {
              if(is_boundary_point_an_extremity[1])
              {
                // the segment is composed of the two boundary points
                return result_type(s);
              }
              else // only boundary_pts[0] is an extremity
              {
                Triangle_3 res_tr = triangle(s.source(), s.target(), vertex(tr, boundary_pts[1]));
                return result_type(std::forward<Triangle_3>(res_tr));
              }
            }
            else // boundary_pts[0] is not an extremity
            {
              if(is_boundary_point_an_extremity[1]) // only boundary_pts[1] is an extremity
              {
                Triangle_3 res_tr = triangle(s.source(), s.target(), vertex(tr, boundary_pts[0]));
                return result_type(std::forward<Triangle_3>(res_tr));
              }
              else // neither boundary points are extremities
              {
                Poly res(4);
                res[0] = vertex(tr, boundary_pts[0]);
                res[1] = vertex(tr, boundary_pts[1]);

                if((res[0] - res[1]) * (s.source() - s.target()) > 0)
                {
                  res[2] = s.target();
                  res[3] = s.source();
                }
                else
                {
                  res[3] = s.target();
                  res[2] = s.source();
                }

                return result_type(std::forward<Poly>(res));
              }
            }
          }

          CGAL_assertion(false);
        }
          break;
        // 2 or 3 segments (and 1 or 2 inside points)
        case 2:
        case 3:
        case 4:
        {
          // @todo do that for a single segment too?
          return build_intersection(tet, tr, points, segments, k);
        }
          break;
        default:
          // can't have more than 4 segments (1 per tet face)
          CGAL_assertion(false);
          break;
      }
    }
      break;

    case 3:
    {
      // the input triangle is entirely contained within the tetrahedron
      return result_type(tr);
    }
      break;
    default:
      CGAL_assertion(false); // never happens (only 3 pts in a tr)
      break;
  }

  return result_type();
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(const typename K::Triangle_3& pl,
             const typename K::Tetrahedron_3& tet,
             const K& k)
{
  return intersection(tet, pl, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
