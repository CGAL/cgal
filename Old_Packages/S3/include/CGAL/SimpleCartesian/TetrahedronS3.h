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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/TetrahedronS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_TETRAHEDRONS3_H
#define CGAL_TETRAHEDRONS3_H

#include <CGAL/SimpleCartesian/PlaneS3.h>
#include <CGAL/solve.h>
#include <vector>
#include <CGAL/predicate_classes_3.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
class TetrahedronS3
{
public:
                 TetrahedronS3() {}
                 TetrahedronS3(const PointS3<FT>& p,
                               const PointS3<FT>& q,
                               const PointS3<FT>& r,
                               const PointS3<FT>& s);

  const PointS3<FT>&    vertex(int i) const;
  const PointS3<FT>&    operator[](int i) const;

  bool           operator==(const TetrahedronS3<FT>& t) const;
  bool           operator!=(const TetrahedronS3<FT>& t) const;
  long           id() const;

  Bbox_3         bbox() const;

  TetrahedronS3<FT> transform(const Aff_transformationS3<FT>& t) const;

  Orientation    orientation() const;
  Oriented_side  oriented_side(const PointS3<FT>& p) const;
  Bounded_side   bounded_side(const PointS3<FT>& p) const;

  bool           has_on_boundary(const PointS3<FT>& p) const;
  bool           has_on_positive_side(const PointS3<FT>& p) const;
  bool           has_on_negative_side(const PointS3<FT>& p) const;
  bool           has_on_bounded_side(const PointS3<FT>& p) const;
  bool           has_on_unbounded_side(const PointS3<FT>& p) const;

  bool           is_degenerate() const;

// private:
  PointS3<FT>    e0;
  PointS3<FT>    e1;
  PointS3<FT>    e2;
  PointS3<FT>    e3;
};


template < class FT >
TetrahedronS3<FT>::TetrahedronS3(const PointS3<FT>& p,
                                 const PointS3<FT>& q,
                                 const PointS3<FT>& r,
                                 const PointS3<FT>& s)
 : e0(p), e1(q), e2(r), e3(s)
{}


template < class FT >
bool
TetrahedronS3<FT>::operator==(const TetrahedronS3<FT>& t) const
{
  if ( orientation() != t.orientation() ) return false;
  std::vector< PointS3<FT> > V1;
  std::vector< PointS3<FT> > V2;
  typename std::vector< PointS3<FT> >::iterator uniq_end1;
  typename std::vector< PointS3<FT> >::iterator uniq_end2;
  int k;
  for ( k=0; k < 4; k++) V1.push_back( vertex(k));
  for ( k=0; k < 4; k++) V2.push_back( t.vertex(k));
  std::sort(V1.begin(), V1.end(), Less_xyz< PointS3<FT> >());
  std::sort(V2.begin(), V2.end(), Less_xyz< PointS3<FT> >());
  uniq_end1 = std::unique( V1.begin(), V1.end());
  uniq_end2 = std::unique( V2.begin(), V2.end());
  V1.erase( uniq_end1, V1.end());
  V2.erase( uniq_end2, V2.end());
  return V1 == V2;
}


template < class FT >
inline
bool
TetrahedronS3<FT>::operator!=(const TetrahedronS3<FT>& t) const
{ return !(*this == t); }


template < class FT >
const PointS3<FT>&
TetrahedronS3<FT>::vertex(int i) const
{
  // modulo 4 is a logical operation, hence cheap
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  switch (i)
    {
     case 0: return e0;
     case 1: return e1;
     case 2: return e2;
    default: return e3;
    }
}


template < class FT >
inline
const PointS3<FT>&
TetrahedronS3<FT>::operator[](int i) const
{ return vertex(i); }

template < class FT >
inline
bool
TetrahedronS3<FT>::has_on_boundary(const PointS3<FT>& p) const
{ return oriented_side(p) == ON_ORIENTED_BOUNDARY; }


template < class FT >
inline
bool
TetrahedronS3<FT>::has_on_positive_side(const PointS3<FT>& p) const
{ return oriented_side(p) == ON_POSITIVE_SIDE; }


template < class FT >
inline
bool
TetrahedronS3<FT>::has_on_negative_side(const PointS3<FT>& p) const
{ return oriented_side(p) == ON_NEGATIVE_SIDE; }


template < class FT >
inline
bool
TetrahedronS3<FT>::has_on_bounded_side(const PointS3<FT>& p) const
{ return oriented_side(p) == ON_BOUNDED_SIDE; }


template < class FT >
inline
bool
TetrahedronS3<FT>::has_on_unbounded_side(const PointS3<FT>& p) const
{ return oriented_side(p) == ON_UNBOUNDED_SIDE; }
        
template < class FT >
Orientation
TetrahedronS3<FT>::orientation() const
{ return CGAL::orientation(vertex(0), vertex(1), vertex(2), vertex(3)); }

template < class FT >
Oriented_side
TetrahedronS3<FT>::oriented_side(const PointS3<FT>& p) const
{
  Orientation o = orientation();
  if (o != ZERO)
    return Oriented_side(o * bounded_side(p));

  CGAL_assertion (!is_degenerate());
  return ON_ORIENTED_BOUNDARY;
}

template < class FT >
Bounded_side
TetrahedronS3<FT>::bounded_side(const PointS3<FT>& p) const
{
  FT alpha, beta, gamma;
  VectorS3<FT> v0 = vertex(1)-vertex(0);
  VectorS3<FT> v1 = vertex(2)-vertex(0);
  VectorS3<FT> v2 = vertex(3)-vertex(0);
  VectorS3<FT> v3 =    p    - vertex(0);
  solve(v0.x(), v0.y(), v0.z(),
        v1.x(), v1.y(), v1.z(),
        v2.x(), v2.y(), v2.z(),
        v3.x(), v3.y(), v3.z(),
        alpha, beta, gamma);

  if (   (alpha < FT(0)) || (beta < FT(0)) || (gamma < FT(0))
      || (alpha + beta + gamma > FT(1)) )
      return ON_UNBOUNDED_SIDE;

  if (   (alpha == FT(0)) || (beta == FT(0)) || (gamma == FT(0))
      || (alpha+beta+gamma == FT(1)) )
    return ON_BOUNDARY;

  return ON_BOUNDED_SIDE;
}

template < class FT >
bool
TetrahedronS3<FT>::is_degenerate() const
{
  PlaneS3<FT> plane(vertex(0), vertex(1), vertex(2));
  return (plane.is_degenerate()) ? true
                                 : plane.has_on_boundary(vertex(3));
}


template < class FT >
inline
Bbox_3
TetrahedronS3<FT>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox()
       + vertex(2).bbox() + vertex(3).bbox();
}


template < class FT >
inline
TetrahedronS3<FT>
TetrahedronS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{
  return TetrahedronS3<FT>(t.transform(vertex(0)),
                           t.transform(vertex(1)),
                           t.transform(vertex(2)),
                           t.transform(vertex(3)));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRONS3
template < class FT >
std::ostream& operator<<(std::ostream& os, const TetrahedronS3<FT>& t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2] << t[3];
    default:
        os << "TetrahedronS3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] ;
        os <<  ", " << t[3] << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRONS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONS3
template < class FT >
std::istream& operator>>(std::istream& is, TetrahedronS3<FT>& t)
{
    PointS3<FT> p, q, r, s;

    is >> p >> q >> r >> s;

    t = TetrahedronS3<FT>(p, q, r, s);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONS3


CGAL_END_NAMESPACE

#endif
