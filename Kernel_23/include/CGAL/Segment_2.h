// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_SEGMENT_2_H
#define CGAL_SEGMENT_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>

namespace CGAL {

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

  typedef Segment_2                           Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Segment_2>::value));

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

  Segment_2(const Point_2 &sp, const Point_2 &ep)
    :  RSegment_2(typename R::Construct_segment_2()(Return_base_tag(), sp,ep)) {}

  typename cpp11::result_of<typename R::Construct_source_2( Segment_2)>::type
  source() const
  { 
    return R_().construct_source_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Construct_target_2( Segment_2)>::type
  target() const
  {
    return R_().construct_target_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Construct_source_2( Segment_2)>::type
  start() const
  {
    return source();
  }

  typename cpp11::result_of<typename R::Construct_target_2( Segment_2)>::type
  end() const
  {
    return target();
  }

  
  typename cpp11::result_of<typename R_::Construct_min_vertex_2(Segment_2)>::type
  min BOOST_PREVENT_MACRO_SUBSTITUTION() const;
  
  typename cpp11::result_of<typename R_::Construct_max_vertex_2( Segment_2)>::type
  max BOOST_PREVENT_MACRO_SUBSTITUTION () const;

  typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2, int)>::type
  vertex(int i) const;

  typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2, int)>::type
  point(int i) const;

  typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2, int)>::type
  operator[](int i) const;

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        has_on(const Point_2 &p) const;  
  bool        collinear_has_on(const Point_2 &p) const;
  FT          squared_length() const;

  bool        is_degenerate() const;

  Bbox_2      bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

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
    return Segment_2(t.transform(source()), t.transform(target()));
  }
};

template < class R_ >
CGAL_KERNEL_INLINE
typename cpp11::result_of<typename R_::Construct_min_vertex_2( Segment_2<R_> )>::type
Segment_2<R_>::min BOOST_PREVENT_MACRO_SUBSTITUTION () const
{
  typename R_::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? source() : target();
}

template < class R_ >
CGAL_KERNEL_INLINE
typename cpp11::result_of<typename R_::Construct_max_vertex_2( Segment_2<R_> )>::type
Segment_2<R_>::max BOOST_PREVENT_MACRO_SUBSTITUTION () const
{
  typename R_::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? target() : source();
}

template < class R_ >
CGAL_KERNEL_INLINE
typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2<R_>, int )>::type
Segment_2<R_>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R_ >
inline
typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2<R_>, int )>::type
Segment_2<R_>::point(int i) const
{
  return vertex(i);
}

template < class R_ >
inline
typename cpp11::result_of<typename R_::Construct_vertex_2( Segment_2<R_>, int )>::type
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



template < class R >
std::ostream &
operator<<(std::ostream &os, const Segment_2<R> &s)
{
    switch(os.iword(IO::mode)) {
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
