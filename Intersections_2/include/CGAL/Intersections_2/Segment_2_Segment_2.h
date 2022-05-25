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
#include <CGAL/Uncertain.h>
#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Cartesian_converter.h>

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
bool
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

    // Let's call s1=[a, B] and s2=[C, D].
    // Then the intersection point I is such that:
    //   I = alpha_s1 × A + (1 - alpha_s1) × B
    //   I = alpha_s2 × C + (1 - alpha_s2) × D
    // or, if `AI` is the notation for the vector `Vector_2(A, I)`:
    //  AI = ( 1 - alpha_s1 ) × AB
    //  CI = ( 1 - alpha_s2 ) × CD
    // If "AB×CD" is the cross product of the two vectors, then
    //
    // Now, let's solve the equations...
    //   AI = AB + BD + DI = AB + BD + (CI - CD)
    // thus
    //   ( 1 - alpha_s1 ) × AB = AB + BD + (1 - alpha_s2)× CD - CD
    //       - alpha_s1   × AB =      BD      - alpha_s2 × CD
    //                                  (that is the equation [1])
    //
    // Let's take the cross product with CD:
    //   -alpha_s1 × ( AB×CD ) = BD × CD
    // and thus:
    //   alpha_s1 = - det(BD, CD) / det(AB, CD) = det(BD, CD) / det(BA, CD)
    //
    // Let's take the cross product of equation [1] with AB:
    //   0 = BD × AB - alpha_s2 × ( CD × AB )
    // and this:
    //   alpha_s2 = det(BD, AB) / det(CD, AB)
    //            = det(BD, BA) / det(CD, BA)
    //            = - det(BD, BA) / det(BA, CD)
    typename K::FT s1_dx = pts[0].x() - pts[1].x(),  // BA
                   s1_dy = pts[0].y() - pts[1].y(),
                   s2_dx = pts[3].x() - pts[2].x(),  // CD
                   s2_dy = pts[3].y() - pts[2].y(),
                   lx    = pts[3].x() - pts[1].x(),  // BD
                   ly    = pts[3].y() - pts[1].y();

    auto det = [](const auto v1, const auto v2) -> typename decltype(v2)::R::FT {
      return v1.x() * v2.y() - v1.y() * v2.x();
    };
#if CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG
#  undef CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG
#  define CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG true
#else
#  undef CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG
#  define CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG false
#endif
    if constexpr (!CGAL_DISABLE_IMPROVEMENT_OF_INTERSECT_SEG_SEG &&
                  std::is_same_v<typename K::FT, double>)
    {
      using Approximate_kernel =
          CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>;
      using Approx_point = CGAL::Point_2<Approximate_kernel>;
      CGAL::Protect_FPU_rounding<true> rounding_mode_protection;
      CGAL::Cartesian_converter<K, Approximate_kernel> convert;
      CGAL::Cartesian_converter<Approximate_kernel, K> convert_back;
      const auto a = convert(pts[0]);
      const auto b = convert(pts[1]);
      const auto c = convert(pts[2]);
      const auto d = convert(pts[3]);
      const auto vector_ba = a - b;
      const auto vector_cd = d - c;
      const auto vector_bd = d - b;
      const auto det_bd_cd = det(vector_bd, vector_cd);
      const auto det_ba_cd = det(vector_ba, vector_cd);
      const auto det_bd_ba = det(vector_bd, vector_ba);
      const auto alpha_s1 = det_bd_cd / det_ba_cd;
      const auto alpha_s2 = - det_bd_ba / det_ba_cd;
      const auto i1 = barycenter(a, alpha_s1, b);
      const auto i2 = barycenter(c, alpha_s2, d);
      const CGAL::Interval_nt_advanced i_x{
          (std::max)(i1.x().inf(), i2.x().inf()),
          (std::min)(i1.x().sup(), i2.x().sup())};
      const CGAL::Interval_nt_advanced i_y{
          (std::max)(i1.y().inf(), i2.y().inf()),
          (std::min)(i1.y().sup(), i2.y().sup())};
      _intersection_point = { std::midpoint(i_x.inf(), i_x.sup()),
                              std::midpoint(i_y.inf(), i_y.sup()) };
      CGAL_assertion(is_valid(_intersection_point.x()));
      CGAL_assertion(is_valid(_intersection_point.y()));
      return _result;
    }
    const auto vector_ba = pts[0] - pts[1];
    const auto vector_cd = pts[3] - pts[2];
    const auto vector_bd = pts[3] - pts[1];
    const auto det_bd_cd = det(vector_bd, vector_cd);
    const auto det_ba_cd = det(vector_ba, vector_cd);
    const auto alpha_s1 = det_bd_cd / det_ba_cd;
    _intersection_point = K().construct_barycenter_2_object()(pts[0], alpha_s1, pts[1]);
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

} //namespace CGAL

#endif
