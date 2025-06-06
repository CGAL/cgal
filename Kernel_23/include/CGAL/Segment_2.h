// Copyright (c) 1999
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_SEGMENT_2_H
#define CGAL_SEGMENT_2_H

#include <CGAL/assertions.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/kernel_config.h>

#include <type_traits>

namespace CGAL {

template <class R_>
class Segment_2 : public R_::Kernel_base::Segment_2
{
  typedef typename R_::Boolean                Boolean;
  typedef typename R_::RT                     RT;
  typedef typename R_::FT                     FT;
  typedef typename R_::Point_2                Point_2;
  typedef typename R_::Direction_2            Direction_2;
  typedef typename R_::Vector_2               Vector_2;
  typedef typename R_::Line_2                 Line_2;
  typedef typename R_::Aff_transformation_2   Aff_transformation_2;
  typedef typename R_::Kernel_base::Segment_2 RSegment_2;

  typedef Segment_2                           Self;
  static_assert(std::is_same<Self, typename R_::Segment_2>::value);

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef RSegment_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_                               R;

  Segment_2() {}

  // conversion from implementation class object to interface class object
  Segment_2(const RSegment_2& s)
    : RSegment_2(s) {}

  Segment_2(RSegment_2&& s)
    : RSegment_2(std::move(s)) {}

  Segment_2(const Point_2 &sp, const Point_2 &ep)
    :  RSegment_2(typename R::Construct_segment_2()(Return_base_tag(), sp,ep)) {}

  decltype(auto)
  source() const
  {
    return R_().construct_source_2_object()(*this);
  }

  decltype(auto)
  target() const
  {
    return R_().construct_target_2_object()(*this);
  }

  decltype(auto)
  start() const
  {
    return source();
  }

  decltype(auto)
  end() const
  {
    return target();
  }

  decltype(auto)
  min BOOST_PREVENT_MACRO_SUBSTITUTION() const {
    typename R_::Less_xy_2 less_xy;
    return less_xy(source(), target()) ? source() : target();
  }

  decltype(auto)
  max BOOST_PREVENT_MACRO_SUBSTITUTION() const {
    typename R_::Less_xy_2 less_xy;
    return less_xy(source(), target()) ? target() : source();
  }

  decltype(auto)
  vertex(int i) const
  { return (i%2 == 0) ? source() : target(); }

  decltype(auto)
  point(int i) const
  { return vertex(i); }

  decltype(auto)
  operator[](int i) const
  { return vertex(i); }

  Boolean is_horizontal() const;
  Boolean is_vertical() const;
  Boolean has_on(const Point_2 &p) const;
  Boolean collinear_has_on(const Point_2 &p) const;
  FT          squared_length() const;

  Boolean is_degenerate() const;

  Bbox_2      bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  Direction_2
  direction() const
  {
    typename R::Construct_vector_2 construct_vector;
    return Direction_2( construct_vector( source(), target()));
  }

  Vector_2
  to_vector() const
  {
    typename R::Construct_vector_2 construct_vector;
    return construct_vector( source(), target());
  }

  Line_2
  supporting_line() const
  {
    typename R::Construct_line_2 construct_line;
    return construct_line(*this);
  }

  Segment_2
  opposite() const
  {
    return R().construct_opposite_segment_2_object()(*this);
  }

  Segment_2
  transform(const Aff_transformation_2 &t) const
  {
    return Segment_2(t.transform(source()), t.transform(target()));
  }
};


template < class R_ >
CGAL_KERNEL_INLINE
typename R_::Boolean
Segment_2<R_>::is_horizontal() const
{
  return R_().equal_y_2_object()(source(), target());
}


template < class R_ >
CGAL_KERNEL_INLINE
typename R_::Boolean
Segment_2<R_>::is_vertical() const
{
  return R_().equal_x_2_object()(source(), target());
}


template < class R_ >
CGAL_KERNEL_INLINE
typename R_::Boolean
Segment_2<R_>::
has_on(const typename R_::Point_2 &p) const
{
  return R_().are_ordered_along_line_2_object()(source(),
                                               p,
                                               target());
}

template < class R_ >
inline
typename R_::Boolean
Segment_2<R_>::
collinear_has_on(const typename R_::Point_2 &p) const
{
  return R_().collinear_has_on_2_object()
               (*this, p);
}

template < class R_ >
CGAL_KERNEL_INLINE
typename Segment_2<R_>::FT
Segment_2<R_>::squared_length() const
{
 return R_().compute_squared_length_2_object()(*this);
}

template < class R_ >
inline
typename R_::Boolean
Segment_2<R_>::is_degenerate() const
{
  return R().is_degenerate_2_object()(*this);
}

template < class R >
std::ostream &
operator<<(std::ostream &os, const Segment_2<R> &s)
{
    switch(IO::get_mode(os)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "Segment_2(" << s.source() <<  ", " << s.target() << ")";
    }
}

template < class R >
std::istream &
operator>>(std::istream &is, Segment_2<R> &s)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
        s = Segment_2<R>(p, q);
    return is;
}

} //namespace CGAL

#endif //  CGAL_SEGMENT_2_H
