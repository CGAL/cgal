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
// file          : include/CGAL/SimpleCartesian/TriangleS3.h
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

#ifndef CGAL_TRIANGLES3_H
#define CGAL_TRIANGLES3_H

#include <CGAL/SimpleCartesian/PlaneS3.h>
#include <CGAL/solve.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
class TriangleS3
{
public:
                 TriangleS3() {}
                 TriangleS3(const PointS3<FT>& p,
                            const PointS3<FT>& q,
                            const PointS3<FT>& r);

  bool           operator==(const TriangleS3<FT>& t) const;
  bool           operator!=(const TriangleS3<FT>& t) const;


  PlaneS3<FT>    supporting_plane() const;

  TriangleS3     transform(const Aff_transformationS3<FT>& t) const;

  bool           has_on(const PointS3<FT>& p) const;
  bool           is_degenerate() const;


  const PointS3<FT>&    vertex(int i) const;
  const PointS3<FT>&    operator[](int i) const;

  Bbox_3         bbox() const;

// private:
  PointS3<FT> e0;
  PointS3<FT> e1;
  PointS3<FT> e2;
};


template < class FT >
TriangleS3<FT>::TriangleS3(const PointS3<FT>& p,
                           const PointS3<FT>& q,
                           const PointS3<FT>& r)
 : e0(p), e1(q), e2(r)
{}


template < class FT >
bool
TriangleS3<FT>::operator==(const TriangleS3<FT>& t) const
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
TriangleS3<FT>::operator!=(const TriangleS3<FT>& t) const
{ return !(*this == t); }


template < class FT >
const PointS3<FT>&
TriangleS3<FT>::vertex(int i) const
{
  if (i<0) i=(i%3)+3;
  else if (i>3) i=i%3;
  return (i==0) ? e0 : (i==1) ? e1 : e2 ;
}


template < class FT >
inline
const PointS3<FT>&
TriangleS3<FT>::operator[](int i) const
{ return vertex(i); }

template < class FT >
inline
PlaneS3<FT>
TriangleS3<FT>::supporting_plane() const
{ return PlaneS3<FT>(vertex(0), vertex(1), vertex(2)); }


template < class FT >
Bbox_3
TriangleS3<FT>::bbox() const
{ return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }


template < class FT >
inline
TriangleS3<FT>
TriangleS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{
  return TriangleS3<FT>(t.transform(vertex(0)),
                        t.transform(vertex(1)),
                        t.transform(vertex(2)));
}
template < class FT >
bool
TriangleS3<FT>::has_on(const PointS3<FT>& p) const
{
  PlaneS3<FT> sp = supporting_plane();
  if ( ! sp.has_on_boundary(p)) return false;
  PointS3<FT> o = vertex(0) + sp.orthogonal_vector();
  FT alpha, beta, gamma;
  VectorS3<FT> v0 = vertex(0)-o;
  VectorS3<FT> v1 = vertex(1)-o;
  VectorS3<FT> v2 = vertex(2)-o;
  VectorS3<FT> v3 =      p  - o;
  solve(v0.x(), v0.y(), v0.z(),
        v1.x(), v1.y(), v1.z(),
        v2.x(), v2.y(), v2.z(),
        v3.x(), v3.y(), v3.z(),
        alpha, beta, gamma);

  return (alpha >= FT(0)) && (beta >= FT(0)) && (gamma >= FT(0))
      && ((alpha+beta+gamma == FT(1)));
}


template < class FT >
bool
TriangleS3<FT>::is_degenerate() const
{
  return collinear(vertex(0),vertex(1),vertex(2));
}


#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLES3
template < class FT >
std::ostream& operator<<(std::ostream& os, const TriangleS3<FT>& t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2];
    default:
        os << "TriangleS3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] <<")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLES3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLES3
template < class FT >
std::istream& operator>>(std::istream& is, TriangleS3<FT>& t)
{
    PointS3<FT> p, q, r;

    is >> p >> q >> r;

    t = TriangleS3<FT>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLES3


CGAL_END_NAMESPACE

#endif
