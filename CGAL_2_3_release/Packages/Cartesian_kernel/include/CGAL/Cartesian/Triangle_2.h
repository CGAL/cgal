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
// file          : include/CGAL/Cartesian/Triangle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_TRIANGLE_2_H
#define CGAL_CARTESIAN_TRIANGLE_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TriangleC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Triangle_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Triangle_handle_2         Triangle_handle_2_;
  typedef typename Triangle_handle_2_::element_type Triangle_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef TriangleC2<R,Cartesian_tag>           Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef TriangleC2<R>                         Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  TriangleC2()
    : Triangle_handle_2_(Triangle_ref_2()) {}

  TriangleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : Triangle_handle_2_(Triangle_ref_2(p, q, r)) {}

  bool           operator==(const Self &s) const;
  bool           operator!=(const Self &s) const;

  Point_2        vertex(int i) const;
  Point_2        operator[](int i) const;

  Self           opposite() const;
  Self           transform(const Aff_transformation_2 &t) const
  {
    return Self(t.transform(vertex(0)),
                t.transform(vertex(1)),
                t.transform(vertex(2)));
  }

  Orientation    orientation() const;
  Oriented_side  oriented_side(const Point_2 &p) const;
  Bounded_side   bounded_side(const Point_2 &p) const;

  bool           has_on_boundary(const Point_2 &p) const;

  bool           has_on_bounded_side(const Point_2 &p) const;
  bool           has_on_unbounded_side(const Point_2 &p) const;

  bool           has_on_positive_side(const Point_2 &p) const;
  bool           has_on_negative_side(const Point_2 &p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;

  FT             area() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleC2<R CGAL_CTAG>::operator==(const TriangleC2<R CGAL_CTAG> &t) const
{
  if (identical(t))
      return true;

  int i;
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
      break;

  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::operator!=(const TriangleC2<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleC2<R CGAL_CTAG>::Point_2
TriangleC2<R CGAL_CTAG>::vertex(int i) const
{
  if (i>2) i = i%3;
  else if (i<0) i = (i%3) + 3;
  return (i==0) ? Ptr()->e0 :
         (i==1) ? Ptr()->e1 :
                  Ptr()->e2;
}

template < class R >
inline
typename TriangleC2<R CGAL_CTAG>::Point_2
TriangleC2<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleC2<R CGAL_CTAG>::FT
TriangleC2<R CGAL_CTAG>::area() const
{
  typename R::Vector_2 v1 = vertex(1)-vertex(0);
  typename R::Vector_2 v2 = vertex(2)-vertex(0);
  return det2x2_by_formula(v1.x(), v1.y(), v2.x(), v2.y())/FT(2);
}

template < class R >
inline
Orientation
TriangleC2<R CGAL_CTAG>::orientation() const
{
  return CGAL::orientation(vertex(0), vertex(1), vertex(2));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
TriangleC2<R CGAL_CTAG>::
bounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  Orientation o1 = CGAL::orientation(vertex(0), vertex(1), p),
              o2 = CGAL::orientation(vertex(1), vertex(2), p),
              o3 = CGAL::orientation(vertex(2), vertex(3), p);

  if (o2 == o1 && o3 == o1)
    return ON_BOUNDED_SIDE;
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_BOUNDARY
     : ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
TriangleC2<R CGAL_CTAG>::
oriented_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  // depends on the orientation of the vertices
  Orientation o1 = CGAL::orientation(vertex(0), vertex(1), p),
              o2 = CGAL::orientation(vertex(1), vertex(2), p),
              o3 = CGAL::orientation(vertex(2), vertex(3), p),
              ot = CGAL::orientation(vertex(0), vertex(1), vertex(2));

  if (o1 == ot && o2 == ot && o3 == ot) // ot cannot be COLLINEAR
    return Oriented_side(ot);
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_ORIENTED_BOUNDARY
     : Oriented_side(-ot);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R CGAL_CTAG>::
has_on_bounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R CGAL_CTAG>::
has_on_unbounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_boundary(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_negative_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_positive_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::is_degenerate() const
{
  return collinear(vertex(0), vertex(1), vertex(2));
}

template < class R >
inline
Bbox_2
TriangleC2<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox();
}

template < class R >
inline
TriangleC2<R CGAL_CTAG>
TriangleC2<R CGAL_CTAG>::
opposite() const
{
  return TriangleC2<R CGAL_CTAG>(vertex(0), vertex(2), vertex(1));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleC2<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleC2(" << t[0] << ", " 
		 << t[1] << ", " << t[2] <<")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2
template < class R >
std::istream &
operator>>(std::istream &is, TriangleC2<R CGAL_CTAG> &t)
{
    TriangleC2<R CGAL_CTAG>::Point_2 p, q, r;

    is >> p >> q >> r;

    if (is)
	t = TriangleC2<R CGAL_CTAG>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRIANGLE_2_H
