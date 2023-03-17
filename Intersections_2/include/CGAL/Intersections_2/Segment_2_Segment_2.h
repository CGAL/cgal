// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_SEGMENT_2_SEGMENT_2_H
#define CGAL_INTERSECTIONS_2_SEGMENT_2_SEGMENT_2_H

#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Intersections_2/Line_2_Line_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace Intersections {

namespace internal {

// indices of lexicographically smallest endpoints
// depending on config parameter of S2S2_inter_info.
// {s1_id0,s1_id1,s2_id0,s2_id1}
constexpr int s2s2_id[8][4] =
{
//  seg1   seg2
  { 0,1,   0,1 },
  { 0,1,   1,0 },
  { 1,0,   0,1 },
  { 1,0,   1,0 },
  { 0,1,   0,1 },
  { 1,0,   0,1 },
  { 0,1,   1,0 },
  { 1,0,   1,0 }
};

// struct used to report the combinaric of the intersection
// of 2 2D segments.
// More information could be gathered if exposed in do_intersect.
// See comments with DI_MORE_INFO_TAG
struct S2S2_inter_info
{
  bool inter = false;
  bool dim = 0;
  std::array<int, 2> pt_ids = {-1,-1};
  // integer in [0,7] indicating segment endpoint ordering for determinism
  // 0: p0 < p1 - p2 < p3 - p0 < p2
  // 1: p0 < p1 - p3 < p2 - p0 < p3
  // 2: p1 < p0 - p2 < p3 - p1 < p2
  // 3: p1 < p0 - p3 < p2 - p1 < p3
  // 4: p2 < p3 - p0 < p1 - p2 < p0
  // 5: p2 < p3 - p1 < p0 - p2 < p1
  // 6: p3 < p2 - p0 < p1 - p3 < p0
  // 7: p3 < p2 - p1 < p0 - p3 < p1
  int config;

  S2S2_inter_info(bool inter, int c=-1)
  : inter(inter)
  , config(c)
  {}

  // intersection is an input endpoint
  S2S2_inter_info(int id)
  : inter(true)
  {
    pt_ids[0]=id;
  }

  // intersection is a segment
  S2S2_inter_info(int id1,int id2)
  : inter(true)
  , dim(1)
  {
    pt_ids[0]=id1;
    pt_ids[1]=id2;
  }
};

template <class K>
inline S2S2_inter_info
do_intersect(const typename K::Segment_2 &seg1, const typename K::Segment_2 &seg2);


// lexicographic order of points p1 < p3 < p2 < p4, with segments (p1,p2) and (p3,p4)
template <class K>
S2S2_inter_info
seg_seg_do_intersect_crossing(
        const typename K::Point_2& p1, const typename K::Point_2& p2,
        const typename K::Point_2& p3, const typename K::Point_2& p4,
        int /* i1 */, int i2, int i3, int /* i4 */,
        const K& k, bool extra_test, int config)
{
    switch (make_certain(k.orientation_2_object()(p1,p2,p3))) {
    case LEFT_TURN:
    {
      switch (k.orientation_2_object()(p3,p4,p2))
      {
        case COLLINEAR:
          return S2S2_inter_info(i2);
        case RIGHT_TURN:
          return S2S2_inter_info(false);
        case LEFT_TURN:
          return S2S2_inter_info(true, config);
        default:
          CGAL_unreachable();
      }
    }
    case RIGHT_TURN:
    {
      switch (k.orientation_2_object()(p3,p4,p2))
      {
        case COLLINEAR:
          return S2S2_inter_info(i2);
        case RIGHT_TURN:
          return S2S2_inter_info(true, config);
        case LEFT_TURN:
          return S2S2_inter_info(false);
        default:
          CGAL_unreachable();
      }
    }
    case COLLINEAR:
      if (extra_test && k.collinear_2_object()(p3,p4,p2))
        return S2S2_inter_info(i3, i2);
      return S2S2_inter_info(i3);
    default:
      CGAL_unreachable();
    }
    CGAL_kernel_assertion(false);
    return S2S2_inter_info(false);
}

// used internally by Arr_segment_traits_2template <class K>
template <class K>
bool
seg_seg_do_intersect_crossing(
        const typename K::Point_2& p1, const typename K::Point_2& p2,
        const typename K::Point_2& p3, const typename K::Point_2& p4,
        const K& k)
{
  return seg_seg_do_intersect_crossing(p1,p2,p3,p4,0,0,0,0,k,false,-1).inter;
}


// lexicographic order of points p1 < p3 < p4 < p2, with segments (p1,p2) and (p3,p4)
template <class K>
S2S2_inter_info
seg_seg_do_intersect_contained(
        const typename K::Point_2& p1, const typename K::Point_2& p2,
        const typename K::Point_2& p3, const typename K::Point_2& p4,
        int /* i1 */, int /* i2 */, int i3, int i4,
        const K& k, bool extra_test, int config)
{
    switch (make_certain(k.orientation_2_object()(p1,p2,p3))) {
    case LEFT_TURN:
    {
      switch (k.orientation_2_object()(p1,p2,p4))
      {
        case COLLINEAR:
          return S2S2_inter_info(i4);
        case RIGHT_TURN:
          return S2S2_inter_info(true, config);
        case LEFT_TURN:
          return S2S2_inter_info(false);
        default:
          CGAL_unreachable();
      }
    }
    case RIGHT_TURN:
    {
      switch (k.orientation_2_object()(p1,p2,p4))
      {
        case COLLINEAR:
          return S2S2_inter_info(i4);
        case RIGHT_TURN:
          return S2S2_inter_info(false);
        case LEFT_TURN:
          return S2S2_inter_info(true, config);
        default:
          CGAL_unreachable();
      }
    }
    case COLLINEAR:
        if (extra_test && k.collinear_2_object()(p3,p4,p2))
          return S2S2_inter_info(i3, i4);
        return S2S2_inter_info(i3);
    default:
      CGAL_unreachable();
    }
    CGAL_kernel_assertion(false);
    return S2S2_inter_info(false);
}

// used internally by Arr_segment_traits_2
template <class K>
bool
seg_seg_do_intersect_contained(
        const typename K::Point_2& p1, const typename K::Point_2& p2,
        const typename K::Point_2& p3, const typename K::Point_2& p4,
        const K& k)
{
  return seg_seg_do_intersect_contained(p1,p2,p3,p4,0,0,0,0,k,false,-1).inter;
}

template <class K>
S2S2_inter_info
do_intersect_with_info(const typename K::Segment_2 &seg1,
                       const typename K::Segment_2 &seg2,
                       const K& k, bool extra_test)
{
    typename K::Less_xy_2 less_xy;

    bool seg1_is_left_to_right = less_xy(seg1.source(),seg1.target());
    bool seg2_is_left_to_right = less_xy(seg2.source(),seg2.target());

    int A1_id = seg1_is_left_to_right ? 0 : 1;
    int A2_id = seg1_is_left_to_right ? 1 : 0;
    int B1_id = seg2_is_left_to_right ? 0 : 1;
    int B2_id = seg2_is_left_to_right ? 1 : 0;

    typename K::Point_2 const & A1 = seg1.point(A1_id);
    typename K::Point_2 const & A2 = seg1.point(A2_id);
    typename K::Point_2 const & B1 = seg2.point(B1_id);
    typename K::Point_2 const & B2 = seg2.point(B2_id);

    typename K::Compare_xy_2 compare_xy;

  // first try to filter using the bbox of the segments
    if (less_xy(A2,B1)
     || less_xy(B2,A1))
        return S2S2_inter_info(false);

    switch(make_certain(compare_xy(A1,B1))) {
    case SMALLER:
        switch(make_certain(compare_xy(A2,B1))) {
        case SMALLER:
            return S2S2_inter_info(false);
        case EQUAL:
            return S2S2_inter_info(A2_id); // DI_MORE_INFO_TAG: A2==B1 but only A2 is reported
        case LARGER:
            switch(make_certain(compare_xy(A2,B2))) {
            case SMALLER:
                return seg_seg_do_intersect_crossing(A1,A2,B1,B2, A1_id,A2_id,B1_id+2,B2_id+2, k, extra_test, (seg1_is_left_to_right ? 0:2) + (seg2_is_left_to_right ? 0:1) );
            case EQUAL:
                // A1 < B1 < B2 = A1
                if (extra_test && k.collinear_2_object()(A1, A2, B1))
                  return S2S2_inter_info(B1_id+2, B2_id+2); // DI_MORE_INFO_TAG: A2==B2 but only B2 is reported
                return S2S2_inter_info(A2_id); // DI_MORE_INFO_TAG: A2==B2 but only A2 is reported
            case LARGER:
                return seg_seg_do_intersect_contained(A1,A2,B1,B2, A1_id,A2_id,B1_id+2,B2_id+2, k, extra_test, (seg1_is_left_to_right ? 0:2) + (seg2_is_left_to_right ? 0:1));
            default:
              CGAL_unreachable();
            }
        default:
          CGAL_unreachable();

        }
    case EQUAL:
        if (extra_test)
        {
          switch(make_certain(compare_xy(A2,B2))) {
          case SMALLER:
            // A1 = B1 < A2 < B2
            if (k.collinear_2_object()(A1,A2,B2))
              return S2S2_inter_info(A1_id, A2_id); // DI_MORE_INFO_TAG: A1==B1 but only A1 is reported
            break;
          case EQUAL:
            // A1 = B1 < A2 = B2
            return S2S2_inter_info(A1_id, A2_id); // DI_MORE_INFO_TAG: A1==B1 and A2==B2 but only A1 and A2 are reported
          case LARGER:
            // A1 = B1 < B2 < A2
            if (k.collinear_2_object()(A1,A2,B2))
              return S2S2_inter_info(B1_id+2, B2_id+2); // DI_MORE_INFO_TAG: A1==B1 but only B1 is reported
          break;
          default:
            CGAL_unreachable();
          }
        }
        return S2S2_inter_info(A1_id); // DI_MORE_INFO_TAG: A1==B1 but only A1 is reported
    case LARGER:
        switch(make_certain(compare_xy(B2,A1))) {
        case SMALLER:
            return S2S2_inter_info(false);
        case EQUAL:
            return S2S2_inter_info(A1_id); // DI_MORE_INFO_TAG: A1==B2 but only A1 is reported
        case LARGER:
            switch(make_certain(compare_xy(B2,A2))) {
            case SMALLER:
                return seg_seg_do_intersect_crossing(B1,B2,A1,A2, B1_id+2,B2_id+2,A1_id,A2_id, k, extra_test, 4 + (seg1_is_left_to_right ? 0:1) + (seg2_is_left_to_right ? 0:2));
            case EQUAL:
                // B1 < A1 < A2 = B2
                if (extra_test && k.collinear_2_object()(B1, A1, B2))
                  return S2S2_inter_info(A1_id, A2_id); // DI_MORE_INFO_TAG: A2==B2 but only A2 is reported
                return S2S2_inter_info(A2_id); // DI_MORE_INFO_TAG: A2==B2 but only A2 is reported
            case LARGER:
                return seg_seg_do_intersect_contained(B1,B2,A1,A2, B1_id+2,B2_id+2,A1_id,A2_id, k, extra_test, 4 + (seg1_is_left_to_right ? 0:1) + (seg2_is_left_to_right ? 0:2));
            default:
              CGAL_unreachable();
            }
        default:
          CGAL_unreachable();
        }
    default:
      CGAL_unreachable();
    }

    CGAL_kernel_assertion(false);
    return S2S2_inter_info(false);
}


template <class K>
typename K::Boolean
do_intersect(const typename K::Segment_2 &seg1,
             const typename K::Segment_2 &seg2,
             const K& k)
{
  return do_intersect_with_info(seg1, seg2, k, false).inter;
}

template <class K>
class Segment_2_Segment_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT, UNKNOWN};
    Segment_2_Segment_2_pair(typename K::Segment_2 const *seg1,
                            typename K::Segment_2 const *seg2)
            : _seg1(seg1), _seg2(seg2) {}

    Intersection_results intersection_type() const;

    typename K::Point_2    intersection_point() const;
    typename K::Segment_2  intersection_segment() const;
protected:
    typename K::Segment_2 const*   _seg1;
    typename K::Segment_2 const *  _seg2;
    mutable Intersection_results       _result = UNKNOWN;
    mutable typename K::Point_2            _intersection_point, _other_point;
};


inline
double s2s2_alpha(double x0, double y0,
                  double x1, double y1,
                  double x2, double y2,
                  double x3, double y3)
{
  const double s1_dx = x0 - x1,
               s1_dy = y0 - y1,
               s2_dx = x3 - x2,
               s2_dy = y3 - y2,
               lx    = x3 - x1,
               ly    = y3 - y1;
  double val = std::fma(lx,s2_dy,-ly*s2_dx)/std::fma(s1_dx,s2_dy,-s1_dy*s2_dx);
  if (val!=val) return 0.5;
  if (val<0) return 0;
  if (val>1) return 1;
  return val;
}

template <class FT>
FT s2s2_alpha(const FT& x0, const FT& y0,
              const FT& x1, const FT& y1,
              const FT& x2, const FT& y2,
              const FT& x3, const FT& y3)
{
  FT s1_dx = x0 - x1,
     s1_dy = y0 - y1,
     s2_dx = x3 - x2,
     s2_dy = y3 - y2,
     lx    = x3 - x1,
     ly    = y3 - y1;
  return (lx*s2_dy-ly*s2_dx)/(s1_dx*s2_dy-s1_dy*s2_dx);
}


template <class K>
typename Segment_2_Segment_2_pair<K>::Intersection_results
Segment_2_Segment_2_pair<K>::intersection_type() const
{
    if (_result!=UNKNOWN)
        return _result;

    S2S2_inter_info inter_info = do_intersect_with_info(*_seg1, *_seg2, K(), true);

    if (!inter_info.inter) {
        _result = NO_INTERSECTION;
        return _result;
    }

    // check if intersection is a segment
    if (inter_info.dim==1)
    {
      _result=SEGMENT;
      _intersection_point = (inter_info.pt_ids[0]>1)
                          ? _seg2->point(inter_info.pt_ids[0]-2)
                          : _seg1->point(inter_info.pt_ids[0]);
      _other_point = inter_info.pt_ids[1]>1
                   ? _seg2->point(inter_info.pt_ids[1]-2)
                   : _seg1->point(inter_info.pt_ids[1]);
      return _result;
    }

    // starting from here we know that the intersection is a point
    _result = POINT;

    // check if intersection is an input endpoint
    if (inter_info.pt_ids[0]>=0)
    {
      _intersection_point = (inter_info.pt_ids[0]>1)
                          ? _seg2->point(inter_info.pt_ids[0]-2)
                          : _seg1->point(inter_info.pt_ids[0]);
      return _result;
    }

    // segments intersect in their interiors
    int c = inter_info.config;
    std::array<typename K::Point_2, 4> pts = c < 4
                                           ? CGAL::make_array( _seg1->point(s2s2_id[c][0]), _seg1->point(s2s2_id[c][1]),
                                                               _seg2->point(s2s2_id[c][2]), _seg2->point(s2s2_id[c][3]) )
                                           : CGAL::make_array( _seg2->point(s2s2_id[c][2]), _seg2->point(s2s2_id[c][3]),
                                                               _seg1->point(s2s2_id[c][0]), _seg1->point(s2s2_id[c][1]) );

    // special case for vertical and horizontal segments
    if (std::is_floating_point<typename K::FT>::value &&
        std::is_same<typename K::Kernel_tag, Cartesian_tag>::value)
    {
      if (pts[0].x()==pts[1].x() && pts[2].y()==pts[3].y())
      {
        _intersection_point = K().construct_point_2_object()(pts[0].x(), pts[2].y());
        return _result;
      }
      if (pts[0].y()==pts[1].y() && pts[2].x()==pts[3].x())
      {
        _intersection_point = K().construct_point_2_object()(pts[2].x(), pts[0].y());
        return _result;
      }
    }

    typename K::FT alpha =  s2s2_alpha(pts[0].x(), pts[0].y(), pts[1].x(), pts[1].y(), pts[2].x(), pts[2].y(), pts[3].x(), pts[3].y());

    _intersection_point = K().construct_barycenter_2_object()(pts[0], alpha, pts[1]);

    return _result;
}


template <class K>
typename K::Point_2
Segment_2_Segment_2_pair<K>::intersection_point() const
{
    if (_result==UNKNOWN)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Segment_2
Segment_2_Segment_2_pair<K>::intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2;
    if (_result==UNKNOWN)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    return Segment_2(_intersection_point, _other_point);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Segment_2>::result_type
intersection(const typename K::Segment_2 &seg1,
             const typename K::Segment_2 &seg2,
             const K&)
{
    typedef Segment_2_Segment_2_pair<K> is_t;
    is_t ispair(&seg1, &seg2);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>(ispair.intersection_segment());
    }
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION_SELF(Segment_2, 2)
CGAL_DO_INTERSECT_FUNCTION_SELF(Segment_2, 2)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_2_SEGMENT_2_SEGMENT_2_H
