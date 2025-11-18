// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb, Laurent Rineau
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_INTERNAL_POLYLINE_H
#define CGAL_MESH_3_INTERNAL_POLYLINE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/iterator.h>

#include <vector>

namespace CGAL {

/// @cond DEVELOPERS
namespace Mesh_3 {
namespace internal {

template <typename Kernel>
class Polyline
{
  typedef typename Kernel::Point_3  Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::FT       FT;

  typedef std::vector<Point_3>      Data;

public:
  typedef typename Data::const_iterator const_iterator;
  typedef std::pair<Point_3, const_iterator> Point_and_location;

  Polyline() {}
  ~Polyline() {}

  /// adds a point at the end of the polyline
  void add_point(const Point_3& p)
  {
    if( points_.empty() || p != end_point() ) {
      points_.push_back(p);
    }
  }

  /// returns the starting point of the polyline
  const Point_3& start_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.front();
  }

  /// returns the ending point of the polyline
  const Point_3& end_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.back();
  }

  /// returns `true` if the polyline is not degenerate
  bool is_valid() const
  {
    return points_.size() > 1;
  }

  /// returns `true` if polyline is a loop
  bool is_loop() const
  {
    return start_point() == end_point();
  }

  const_iterator next(const_iterator it, Orientation orientation) const {
    if(orientation == POSITIVE) {
      CGAL_assertion(it != (points_.end() - 1));
      if(it == (points_.end() - 2)) {
        CGAL_assertion(is_loop());
        it = points_.begin();
      } else {
        ++it;
      }
    } else {
      CGAL_assertion(orientation == NEGATIVE);
      CGAL_assertion(it != points_.begin());
      if(it == (points_.begin() + 1)) {
        CGAL_assertion(is_loop());
        it = points_.end() - 1;
      } else {
        --it;
      }
    }
    return it;
  }

  bool is_curve_segment_covered(CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1, const FT sq_r2,
                                const const_iterator cc1_it,
                                const const_iterator cc2_it) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);
    typename Kernel::Has_on_bounded_side_3 cover_pred =
      Kernel().has_on_bounded_side_3_object();

    typedef typename Kernel::Sphere_3 Sphere_3;
    const Sphere_3 s1(c1, sq_r1);
    const Sphere_3 s2(c2, sq_r2);

    const_iterator c1_it = cc1_it;
    const_iterator c2_it = cc2_it;

    if(orientation == CGAL::NEGATIVE) {
      ++c1_it;
      ++c2_it;
      CGAL_assertion(c1_it != points_.end());
      CGAL_assertion(c2_it != points_.end());
    }

    if(c1_it == c2_it) return cover_pred(s1, s2, c1, c2);
    const_iterator next_it = this->next(c1_it, orientation);

    if(!cover_pred(s1, s2, c1, *next_it)) return false;

    for(const_iterator it = next_it; it != c2_it; /* in body */) {
      next_it = this->next(it, orientation);
      if(!cover_pred(s1, s2, *it, *next_it)) return false;
      it = next_it;
    } // end loop ]c1_it, c2_it[

    return cover_pred(s1, s2, *c2_it, c2);
  }

  bool is_curve_segment_covered(CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1,    const FT sq_r2) const
  {
    return is_curve_segment_covered(orientation,
                                    c1, c2, sq_r1, sq_r2, locate(c1), locate(c2));
  }

  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          CGAL::Orientation orientation) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);
    const_iterator p_it = locate(p);
    const_iterator q_it = locate(q);
    return curve_segment_length(p, q, orientation, p_it, q_it);
  }

  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          CGAL::Orientation orientation,
                          const_iterator p_it,
                          const_iterator q_it) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);
    CGAL_assertion(p_it == locate(p));
    CGAL_assertion(q_it == locate(q));

    if(p_it == q_it) {
      const CGAL::Comparison_result cmp = compare_distance(*p_it,p,q);
      if( (cmp != LARGER  && orientation == POSITIVE) ||
          (cmp != SMALLER && orientation == NEGATIVE) )
      {
        // If the orientation of `p` and `q` on the segment is compatible
        // with `orientation`, then return the distance between the two
        // points.
        return distance(p, q);
      }
    }

    if(orientation == CGAL::NEGATIVE) {
      ++p_it;
      ++q_it;
      CGAL_assertion(p_it != points_.end());
      CGAL_assertion(q_it != points_.end());
    }

    const_iterator next_it = this->next(p_it, orientation);
    FT result = distance(p, *next_it);
    for(const_iterator it = next_it; it != q_it; /* in body */) {
      next_it = this->next(it, orientation);
      result += distance(*it, *next_it);
      it = next_it;
    } // end loop ]p_it, q_it[
    result += distance(*q_it, q);
    return result;
  }


  /// returns the angle at the first point.
  /// \pre The polyline must be a loop.
  Angle angle_at_first_point() const {
    CGAL_precondition(is_loop());
    const Point_3& first = points_.front();
    const Point_3& next_p = points_[1];
    const Point_3& prev = points_[points_.size() - 2];
    return angle(prev, first, next_p);
  }

  /// returns the length of the polyline
  FT length() const
  {
    if(length_ < 0.)
    {
      FT result(0);
      const_iterator it = points_.begin();
      const_iterator previous = it++;

      for(const_iterator end = points_.end(); it != end; ++it, ++previous) {
        result += distance(*previous, *it);
      }
      length_ = result;
    }
    return length_;
  }

  /// returns the signed geodesic distance between `p` and `q`.
  FT signed_geodesic_distance(const Point_3& p, const Point_3& q) const
  {
    // Locate p & q on polyline
    const_iterator pit = locate(p);
    const_iterator qit = locate(q,false);

    return signed_geodesic_distance(p, q, pit, qit);
  }

  FT signed_geodesic_distance(const Point_3& p, const Point_3& q,
                              const_iterator pit, const_iterator qit) const
  {
    CGAL_assertion(pit == locate(p));
    CGAL_assertion(qit == locate(q, false));

    // If p and q are in the same segment of the polyline
    if ( pit == qit )
    {
      FT result = distance(p,q);

      // Find the closest point to *pit
      if ( compare_distance(*pit,p,q) != CGAL::LARGER )
      { return result; }
      else
      { return -result; }
    }
    if(is_loop())
    {
      FT positive_distance, negative_distance;
      if(pit <= qit)
      {
        positive_distance = curve_segment_length(p, q, CGAL::POSITIVE, pit, qit);
        negative_distance = length() - positive_distance;
      }
      else
      {
        negative_distance = curve_segment_length(q, p, CGAL::POSITIVE, qit, pit);
        positive_distance = length() - negative_distance;
      }
      return (positive_distance < negative_distance) ? positive_distance : (-negative_distance);
    }
    else
    {
      return (pit <= qit)
        ?     curve_segment_length(p, q, CGAL::POSITIVE, pit, qit)
        : ( - curve_segment_length(p, q, CGAL::NEGATIVE, pit, qit) );
    }
  }

  const_iterator previous_segment_source(const_iterator it) const
  {
    CGAL_assertion(it != points_.end());
    if(it == first_segment_source())
    {
      CGAL_assertion(is_loop());
      it = last_segment_source();
    }
    else
    {
      --it;
    }
    return it;
  }

  const_iterator next_segment_source(const_iterator it) const
  {
    CGAL_assertion(it != points_.end());
    if(it == last_segment_source())
    {
      if(is_loop())
        return first_segment_source();
      else
        return points_.end();
    }
    else
    {
      ++it;
      return it;
    }
  }

  /// returns a point at geodesic distance `distance` from p along the
  /// polyline. The polyline is oriented from starting point to end point.
  /// The distance could be negative.
  Point_3 point_at(const Point_3& start_pt, FT distance) const
  {
    return point_at(start_pt, distance, locate(start_pt)).first;
  }

  /// returns a point at geodesic distance `distance` from `start_pt` along the
  /// polyline. The polyline is oriented from starting point to end point.
  /// The distance could be negative.
  Point_and_location point_at(const Point_3& start_pt,
                              FT distance,
                              const_iterator start_it) const
  {
    CGAL_assertion(start_it == locate(start_pt));

    const Point_3& start_it_pt = *start_it;
    const_iterator start_it_locate_pt
      = (start_it == points_.begin()) ? start_it : std::prev(start_it);

    distance += curve_segment_length(start_it_pt, start_pt, CGAL::POSITIVE,
                                     start_it_locate_pt, start_it);

    // If polyline is a loop, ensure that distance is given from start_it
    if(is_loop())
    {
      if(distance < FT(0))         { distance += length(); }
      else if(distance > length()) { distance -= length(); }
    }
    else if(distance < FT(0)) // If polyline is not a loop and distance is < 0, go backward
    {
      Point_3 new_start = start_pt;
      while(distance < FT(0))
      {
        start_it = previous_segment_source(start_it);
        distance += this->distance(new_start, *start_it);
        new_start = *start_it;
      }
    }

    CGAL_assertion(distance >= FT(0));
    CGAL_assertion(distance <= length());

    // Initialize iterators
    const_iterator pit = start_it; // start at start_it, and walk forward
    const_iterator previous = pit++;

    // Iterate to find which segment contains the point we want to construct
    FT segment_length = this->distance(*previous, *pit);
    while(distance > segment_length)
    {
      distance -= segment_length;

      // Increment iterators and update length
      previous = next_segment_source(previous);
      pit = next_segment_source(pit);

      if(pit == points_.end())
      {
        CGAL_assertion(distance < this->distance(*previous, end_point()));
        break; // return {*previous, previous}
      }
      else
        segment_length = this->distance(*previous, *pit);
    }

    // return point at distance from current segment source
    using Vector_3 = typename Kernel::Vector_3;
    auto vector = Kernel().construct_vector_3_object();
    Vector_3 v = (pit != points_.end()) ? vector(*previous, *pit)
                                        : vector(*previous, end_point());

    return {(*previous) + (distance / CGAL::sqrt(v.squared_length())) * v,
             previous};
  }

  const_iterator locate_corner(const Point_3& p) const
  {
    const_iterator res = points_.end();
    if(p == start_point())
      res = points_.begin();
    else if(p == end_point())
      res = last_segment_source();

    CGAL_assertion(res == locate(p));
    CGAL_assertion(res != points_.end());
    return res;
  }

  const_iterator locate_point(const Point_3& p) const
  {
    return locate(p);
  }

  bool are_ordered_along(const Point_3& p, const Point_3& q,
                         const_iterator pit, const_iterator qit) const
  {
    CGAL_precondition(!is_loop());

    // Locate p & q on polyline
    CGAL_assertion(pit == locate(p));
    CGAL_assertion(qit == locate(q, true));

    // Points are not located on the same segment
    if ( pit != qit ) { return (pit <= qit); }

    // pit == qit, then we have to sort p&q along (pit,pit+1)
    return ( compare_distance(*pit,p,q) != CGAL::LARGER );
  }

  bool are_ordered_along(const Point_3& p, const Point_3& q) const
  {
    return are_ordered_along(p, q, locate(p), locate(q, true));
  }

private:
  const_iterator first_segment_source() const
  {
    CGAL_precondition(is_valid());
    return points_.begin();
  }

  const_iterator last_segment_source() const
  {
    CGAL_precondition(is_valid());
    return (points_.end() - 2);
  }

  /// returns an iterator on the starting point of the segment of the
  /// polyline which contains p
  /// if end_point_first is true, then --end is returned instead of begin
  /// if p is the starting point of a loop.
  const_iterator locate(const Point_3& p, bool end_point_first=false) const
  {
    CGAL_precondition(is_valid());

    // First look if p is one of the points of the polyline
    const_iterator result = std::find(points_.begin(), points_.end(), p);
    if ( result != points_.end() )
    {
      if ( result != points_.begin() )
      { return --result; }
      else
      {
        // Treat loops
        if ( end_point_first && p == end_point() )
        { return last_segment_source(); }
        else
        { return result; }
      }
    }

    CGAL_assertion(result == points_.end());

    // Get result by projecting p on the polyline
    const_iterator it = points_.begin();
    const_iterator previous = it;
    Segment_3 nearest_segment;
    const_iterator nearest_vertex = it;
    result = nearest_vertex;
    bool nearest_is_a_segment = false;

    while ( ++it != points_.end() )
    {
      Segment_3 seg (*previous, *it);

      if(nearest_is_a_segment)
      {
        if(compare_distance(p, *it, nearest_segment) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          nearest_is_a_segment = false;
          result = it;
          if (possibly(angle(*previous, *it, p) == CGAL::ACUTE) &&
              compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
          {
            nearest_segment = seg;
            nearest_is_a_segment = true;
            result = previous;
          }
        }
        else if(compare_distance(p, seg, nearest_segment) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          result = previous;
        }
      }
      else {
        if(compare_distance(p, *it, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          result = it;
        }
        if ((nearest_vertex != it ||
             possibly(angle(*previous, *it, p) == CGAL::ACUTE)) &&
            compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          nearest_is_a_segment = true;
          result = previous;
        }
      }
      previous = it;
    } // end the while loop on the vertices of the polyline

    if(result == points_.begin()) {
      return (end_point_first && !nearest_is_a_segment) ? last_segment_source() : points_.begin();
    } else {
      return result;
    }
  }


  FT distance(const Point_3& p, const Point_3& q) const
  {
    typename Kernel::Compute_squared_distance_3 sq_distance =
      Kernel().compute_squared_distance_3_object();
    return CGAL::sqrt(sq_distance(p, q));
  }

  Angle angle(const Point_3& p,
              const Point_3& angle_vertex_point,
              const Point_3& q) const
  {
    typename Kernel::Angle_3 compute_angle =  Kernel().angle_3_object();
    return compute_angle(p,angle_vertex_point,q);
  }

  template <typename T1, typename T2>
  CGAL::Sign compare_distance(const Point_3& p,
                              const T1& obj1,
                              const T2& obj2) const
  {
    typename Kernel::Compare_distance_3 compare_distance =
      Kernel().compare_distance_3_object();
    return compare_distance(p,obj1,obj2);
  }

public:
  Data points_;

private:
  mutable FT length_ = -1.;

}; // end class Polyline


} // end namespace internal
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_INTERNAL_POLYLINE_H
