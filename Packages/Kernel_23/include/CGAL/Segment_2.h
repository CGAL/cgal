// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_SEGMENT_2_H
#define CGAL_SEGMENT_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_2 : public R_::Kernel_base::Segment_2
{
  typedef typename R_::RT                     RT;
  typedef typename R_::FT                     FT;
  typedef typename R_::Point_2                Point_2;
  typedef typename R_::Direction_2            Direction_2;
  typedef typename R_::Vector_2               Vector_2;
  typedef typename R_::Line_2                 Line_2;
  typedef typename R_::Aff_transformation_2   Aff_transformation_2;
  typedef typename R_::Kernel_base::Segment_2 RSegment_2;

public:

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

  Segment_2() 
    : RSegment_2(typename R::Construct_segment_2()().rep())
   {}

  Segment_2(const Point_2 &sp, const Point_2 &ep)
    :  RSegment_2(typename R::Construct_segment_2()(sp,ep).rep()) {}

  // conversion from implementation class object to interface class object
  Segment_2(const RSegment_2& s)
    : RSegment_2(s) {}

  const Point_2&     source() const
  {
    return R_().construct_source_2_object()(*this);
  }

  const Point_2&     target() const
  {
    return R_().construct_target_2_object()(*this);
  }

  const Point_2 &    start() const;
  const Point_2 &    end() const;

  const Point_2 &   min() const;
  const Point_2 &   max() const;
  const Point_2 &   vertex(int i) const;
  const Point_2 &   point(int i) const;
  const Point_2 &   operator[](int i) const;

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        has_on(const Point_2 &p) const;  
  bool        collinear_has_on(const Point_2 &p) const;
  FT          squared_length() const;

  bool        is_degenerate() const;

  Bbox_2      bbox() const;

  bool
  operator==(const Segment_2 &s) const
  {
    return R().equal_2_object()(*this, s);
  }


  bool
  operator!=(const Segment_2 &s) const
  {
    return !(*this == s);
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
    return rep().transform(t);
  }
};

template < class R_ >
inline
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::start() const
{
  return source();
}

template < class R_ >
inline
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::end() const
{
  return target();
}

template < class R_ >
CGAL_KERNEL_INLINE
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::min() const
{
  typename R_::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? source() : target();
}

template < class R_ >
CGAL_KERNEL_INLINE
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::max() const
{
  typename R_::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? target() : source();
}

template < class R_ >
CGAL_KERNEL_INLINE
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R_ >
inline
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::point(int i) const
{
  return vertex(i);
}

template < class R_ >
inline
const typename Segment_2<R_>::Point_2 &
Segment_2<R_>::operator[](int i) const
{
  return vertex(i);
}

template < class R_ >
CGAL_KERNEL_INLINE
bool
Segment_2<R_>::is_horizontal() const
{
  return R_().equal_y_2_object()(source(), target());
}


template < class R_ >
CGAL_KERNEL_INLINE
bool
Segment_2<R_>::is_vertical() const
{
  return R_().equal_x_2_object()(source(), target());
}


template < class R_ >
CGAL_KERNEL_INLINE
bool
Segment_2<R_>::
has_on(const typename R_::Point_2 &p) const
{
  return R_().are_ordered_along_line_2_object()(source(), 
					       p, 
					       target());
}


template < class R_ >
inline
bool
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
bool
Segment_2<R_>::is_degenerate() const
{
  return R().is_degenerate_2_object()(*this);
}


template < class R_ >
CGAL_KERNEL_INLINE
Bbox_2
Segment_2<R_>::bbox() const
{
  return R().construct_bbox_2_object()(*this);
}



#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_2
template < class R>
std::ostream &
operator<<(std::ostream &os, const Segment_2<R> &s)
{
  return os << s.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2
template < class R>
std::istream &
operator>>(std::istream &is, Segment_2<R> &s)
{
  return is >> s.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2

CGAL_END_NAMESPACE

#endif //  CGAL_SEGMENT_2_H
