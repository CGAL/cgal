// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/TriangleS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_TRIANGLES2_H
#define CGAL_TRIANGLES2_H

#include <CGAL/SimpleCartesian/predicates_on_pointsS2.h>

CGAL_BEGIN_NAMESPACE



template <class FT>
class TriangleS2
{
public:
                 TriangleS2() {}
                 TriangleS2(const PointS2<FT>& p,
                             const PointS2<FT>& q,
                             const PointS2<FT>& r);

  bool           operator==(const TriangleS2<FT>& s) const;
  bool           operator!=(const TriangleS2<FT>& s) const;

  const PointS2<FT>&    vertex(int i) const;
  const PointS2<FT>&    operator[](int i) const;

  TriangleS2<FT> transform(const Aff_transformationS2<FT>& t) const;

  Orientation    orientation() const;
  Oriented_side  oriented_side(const PointS2<FT>& p) const;
  Bounded_side   bounded_side(const PointS2<FT>& p) const;

  bool           has_on_boundary(const PointS2<FT>& p) const;

  bool           has_on_bounded_side(const PointS2<FT>& p) const;
  bool           has_on_unbounded_side(const PointS2<FT>& p) const;

  bool           has_on_positive_side(const PointS2<FT>& p) const;
  bool           has_on_negative_side(const PointS2<FT>& p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;

// private:
  PointS2<FT> e0;
  PointS2<FT> e1;
  PointS2<FT> e2;
};


template < class FT >
CGAL_KERNEL_CTOR_INLINE
TriangleS2<FT>::TriangleS2(const PointS2<FT>& p,
                             const PointS2<FT>& q,
                             const PointS2<FT>& r)
 : e0(p), e1(q), e2(r) {}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool 
TriangleS2<FT>::operator==(const TriangleS2<FT>& t) const
{
  int i;
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
      break;

  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class FT >
inline
bool 
TriangleS2<FT>::operator!=(const TriangleS2<FT>& t) const
{ return !(*this == t); }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
const PointS2<FT>&
TriangleS2<FT>::vertex(int i) const
{
  if (i>2) i = i%3;
  else if (i<0) i = (i%3) + 3;
  return (i==0) ? e0 :
         (i==1) ? e1 :
                  e2 ;
}

template < class FT >
inline
const PointS2<FT>&
TriangleS2<FT>::operator[](int i) const
{ return vertex(i); }

template < class FT >
inline
Orientation
TriangleS2<FT>::orientation() const
{ return CGAL::orientation(e0,e1,e2); }

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
TriangleS2<FT>::bounded_side(const PointS2<FT>& p) const
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


template < class FT >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
TriangleS2<FT>::oriented_side(const PointS2<FT>& p) const
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
     : Oriented_side(opposite(ot));
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleS2<FT>::has_on_bounded_side(const PointS2<FT>& p) const
{ return bounded_side(p) == ON_BOUNDED_SIDE; }

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleS2<FT>::has_on_unbounded_side(const PointS2<FT>& p) const
{ return bounded_side(p) == ON_UNBOUNDED_SIDE; }

template < class FT >
inline
bool 
TriangleS2<FT>::has_on_boundary(const PointS2<FT>& p) const
{ return bounded_side(p) == ON_BOUNDARY; }

template < class FT >
inline
bool 
TriangleS2<FT>::has_on_negative_side(const PointS2<FT>& p) const
{ return oriented_side(p) == ON_NEGATIVE_SIDE; }

template < class FT >
inline
bool 
TriangleS2<FT>::has_on_positive_side(const PointS2<FT>& p) const
{ return oriented_side(p) == ON_POSITIVE_SIDE; }

template < class FT >
inline
bool 
TriangleS2<FT>::is_degenerate() const
{ return collinear(vertex(0), vertex(1), vertex(2)); }

template < class FT >
inline
Bbox_2 
TriangleS2<FT>::bbox() const
{ return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }


template < class FT >
inline
TriangleS2<FT>
TriangleS2<FT>::transform(const Aff_transformationS2<FT>& t) const
{
  return TriangleS2<FT>(t.transform(vertex(0)),
                        t.transform(vertex(1)),
                        t.transform(vertex(2)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLES2
template < class FT >
std::ostream& operator<<(std::ostream &os, const TriangleS2<FT> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleS2(" << t[0] << ", " << t[1] << ", " << t[2] <<")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLES2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLES2
template < class FT >
std::istream& operator>>(std::istream &is, TriangleS2<FT> &t)
{
    PointS2<FT> p, q, r;

    is >> p >> q >> r;

    t = TriangleS2<FT>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLES2



CGAL_END_NAMESPACE

#endif
