// Copyright (c) 2009  GeometryFactory (France), INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//
#ifndef CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Line_3_intersection.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type
intersection_collinear_segments(const typename K::Segment_3& s1,
                                const typename K::Segment_3& s2,
                                const K& k)
{
  CGAL_precondition(!s1.is_degenerate() && !s2.is_degenerate());

  const typename K::Point_3& p=s1[0], q=s1[1], r=s2[0], s=s2[1];

  typename K::Collinear_are_ordered_along_line_3 cln_order=k.collinear_are_ordered_along_line_3_object();

  if(cln_order(p,r,q))
  {
    if(cln_order(p,s,q))
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s2);

    if(cln_order(r,p,s))
    {
      if(r != p)
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(r,p) );

      if(cln_order(r,q,s))
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1);

      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    }

    return (r != q) ? intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(r,q) )
                    : intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(q);
  }

  if(cln_order(p,s,q))
  {
    if(cln_order(r,p,s))
    {
      if(s != p)
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(s,p) );
      if(cln_order(r,q,s))
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1);

      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    }

    return (s != q) ? intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(s,q) )
                    : intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(q);
  }

  if(cln_order(r,p,s))
    return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1);

  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
}

template<class K>
struct L_p_visitor
{
  typedef typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type result_type;

  L_p_visitor(const typename K::Segment_3& s1,
              const typename K::Segment_3& s2)
    : s1(s1), s2(s2)
  { }

  const typename K::Segment_3& s1;
  const typename K::Segment_3& s2;

  result_type operator()(const typename K::Point_3& p) const
  {
    typename K::Collinear_are_ordered_along_line_3 cln_order = K().collinear_are_ordered_along_line_3_object();

    if(cln_order(s1[0],p,s1[1]) && cln_order(s2[0],p,s2[1]))
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    else
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
  }

  result_type operator()(const typename K::Line_3&) const
  {
    return intersection_collinear_segments(s1,s2,K());
  }
};

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type
intersection(const typename K::Segment_3& s1,
             const typename K::Segment_3& s2,
             const K& k)
{
  CGAL_precondition(!s1.is_degenerate() && !s2.is_degenerate());

  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type
      v = internal::intersection(s1.supporting_line(), s2.supporting_line(), k);

  if(v)
    return std::visit(L_p_visitor<K>(s1, s2) , *v);

  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_SEGMENT_3_INTERSECTION_H
