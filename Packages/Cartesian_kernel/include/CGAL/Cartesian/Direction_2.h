// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Direction_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_2_H
#define CGAL_CARTESIAN_DIRECTION_2_H

#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class DirectionC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Direction_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Direction_handle_2        Direction_handle_2_;
  typedef typename Direction_handle_2_::element_type Direction_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef DirectionC2<R,Cartesian_tag>          Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef DirectionC2<R>                        Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  DirectionC2()
    : Direction_handle_2_(Direction_ref_2()) {}

  DirectionC2(const Vector_2 &v)
    : Direction_handle_2_(v) {}

  DirectionC2(const Line_2 &l)
    : Direction_handle_2_(l.direction()) {}

  DirectionC2(const Ray_2 &r)
    : Direction_handle_2_(r.direction()) {}

  DirectionC2(const Segment_2 &s)
    : Direction_handle_2_(s.direction()) {}

  DirectionC2(const FT &x, const FT &y)
    : Direction_handle_2_(Direction_ref_2(x, y)) {}

  bool operator==(const Self &d) const;
  bool operator!=(const Self &d) const;
  bool operator>=(const Self &d) const;
  bool operator<=(const Self &d) const;
  bool operator>(const Self &d) const;
  bool operator<(const Self &d) const;
  bool counterclockwise_in_between( const Self &d1, const Self &d2) const;
  
  Vector_2 to_vector() const;

  Self perpendicular(const Orientation &o) const;
  Self transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

  Self operator-() const;

  FT delta(int i) const;
  FT dx() const
  {
      return Ptr()->e0;
  }
  FT dy() const
  {
      return Ptr()->e1;
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
bool
DirectionC2<R CGAL_CTAG>::operator==(const DirectionC2<R CGAL_CTAG> &d) const
{
  if (identical(d))
      return true;
  return equal_direction(*this, d);
}

template < class R >
inline
bool
DirectionC2<R CGAL_CTAG>::operator!=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return !( *this == d );
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator<(const DirectionC2<R CGAL_CTAG> &d) const
{
  return compare_angle_with_x_axis(*this, d) == SMALLER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator>(const DirectionC2<R CGAL_CTAG> &d) const
{
  return d < *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator>=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return compare_angle_with_x_axis(*this, d) != SMALLER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator<=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return compare_angle_with_x_axis(*this, d) != LARGER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::
counterclockwise_in_between(const DirectionC2<R CGAL_CTAG> &d1,
                            const DirectionC2<R CGAL_CTAG> &d2) const
// returns true, iff \ccVar\ is not equal to \ccc{d1}, and 
// while rotating counterclockwise starting at \ccc{d1}, 
// \ccVar\ is reached strictly before \ccc{d2} is reached.
// Note that true is returned if \ccc{d1} == \ccc{d2}, unless
//  also \ccVar\ == \ccc{d1}.
{
  if ( d1 < *this)
    return ( *this < d2 )||( d2 <= d1 );
  else
    return ( *this < d2 )&&( d2 <= d1 );
}

template < class R >
inline
typename DirectionC2<R CGAL_CTAG>::Vector_2
DirectionC2<R CGAL_CTAG>::to_vector() const
{
  return Vector_2(*this);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
DirectionC2<R CGAL_CTAG>
DirectionC2<R CGAL_CTAG>::perpendicular(const Orientation &o) const
{ // FIXME : construction
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
    return DirectionC2<R CGAL_CTAG>(-dy(), dx());
  else
    return DirectionC2<R CGAL_CTAG>(dy(), -dx());
}

template < class R >
inline
DirectionC2<R CGAL_CTAG>
DirectionC2<R CGAL_CTAG>::operator-() const
{ // FIXME : construction
  return DirectionC2<R CGAL_CTAG>(-dx(), -dy());
}

template < class R >
CGAL_KERNEL_INLINE
typename DirectionC2<R CGAL_CTAG>::FT
DirectionC2<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  return (i==0) ? dx() : dy();
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC2
template < class R >
std::ostream&
operator<<(std::ostream &os, const DirectionC2<R CGAL_CTAG> &d)
{
    typename R::Vector_2 v = d.to_vector();
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "DirectionC2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2
template < class R >
std::istream&
operator>>(std::istream &is, DirectionC2<R CGAL_CTAG> &p)
{
    typename DirectionC2<R CGAL_CTAG>::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << std::endl << "Stream must be in ascii or binary mode"
	          << std::endl;
        break;
    }
    if (is)
	p = DirectionC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DIRECTION_2_H
