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
// file          : include/CGAL/SimpleCartesian/PlaneS3.h
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

#ifndef CGAL_PLANES3_H
#define CGAL_PLANES3_H

#include <CGAL/SimpleCartesian/PointS2.h>
#include <CGAL/solve.h>
#include <CGAL/SimpleCartesian/basic_constructionsS3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class PlaneS3
{
public:
                   PlaneS3() {}
                   PlaneS3(const PointS3<FT>& p,
                           const PointS3<FT>& q,
                           const PointS3<FT>& r);
                   PlaneS3(const PointS3<FT>& p,
                           const DirectionS3<FT>& d);
                   PlaneS3(const PointS3<FT>& p,
                           const VectorS3<FT>& v);
                   PlaneS3(const FT& a, const FT& b,
                           const FT& c, const FT& d);
                   PlaneS3(const LineS3<FT>& l,
                           const PointS3<FT>& p);
                   PlaneS3(const SegmentS3<FT>& s,
                           const PointS3<FT>& p);
                   PlaneS3(RayS3<FT>& r,
                           const PointS3<FT>& p);

  bool             operator==(const PlaneS3<FT>& p) const;
  bool             operator!=(const PlaneS3<FT>& p) const;

  const FT&        a() const;
  const FT&        b() const;
  const FT&        c() const;
  const FT&        d() const;

  LineS3<FT>       perpendicular_line(const PointS3<FT>& p) const;
  PlaneS3          opposite() const;

  PointS3<FT>      point() const;
  PointS3<FT>      projection(const PointS3<FT>& p) const;
  VectorS3<FT>     orthogonal_vector() const;
  DirectionS3<FT>  orthogonal_direction() const;
  VectorS3<FT>     base1() const;
  VectorS3<FT>     base2() const;

  PointS3<FT>      to_plane_basis(const PointS3<FT>& p) const;

  PointS2<FT>      to_2d(const PointS3<FT>& p) const;
  PointS3<FT>      to_3d(const PointS2<FT>& p) const;

  PlaneS3          transform(const Aff_transformationS3<FT>& t) const;


  Oriented_side    oriented_side(const PointS3<FT>& p) const;
  bool             has_on_boundary(const PointS3<FT>& p) const;
  bool             has_on_boundary(const LineS3<FT>& p) const;
  bool             has_on_positive_side(const PointS3<FT>& l) const;
  bool             has_on_negative_side(const PointS3<FT>& l) const;

  bool             is_degenerate() const;

// private:
  void              new_rep(const PointS3<FT>& p,
                            const PointS3<FT>& q,
                            const PointS3<FT>& r);
  void              new_rep(const FT& a, const FT& b,
                            const FT& c, const FT& d);

  FT    e0;
  FT    e1;
  FT    e2;
  FT    e3;
};


template < class FT >
inline
void
PlaneS3<FT>::new_rep(const FT& a, const FT& b, const FT& c, const FT& d)
{
  e0 = a;
  e1 = b;
  e2 = c;
  e3 = d;
}

template < class FT >
inline
void
PlaneS3<FT>::new_rep(const PointS3<FT>& p,
                     const PointS3<FT>& q,
                     const PointS3<FT>& r)
{
  FT rpx = p.x()-r.x();
  FT rpy = p.y()-r.y();
  FT rpz = p.z()-r.z();
  FT rqx = q.x()-r.x();
  FT rqy = q.y()-r.y();
  FT rqz = q.z()-r.z();
  // Cross product rp * rq.
  e0 = rpy*rqz - rqy*rpz;
  e1 = rpz*rqx - rqz*rpx;
  e2 = rpx*rqy - rqx*rpy;
  e3 = - e0*r.x() - e1*r.y() - e2*r.z();
}


CGAL_END_NAMESPACE

#include <CGAL/SimpleCartesian/LineS3.h>

CGAL_BEGIN_NAMESPACE


template < class FT >
PlaneS3<FT>::PlaneS3(const PointS3<FT>& p,
                     const PointS3<FT>& q,
                     const PointS3<FT>& r)
{ new_rep(p, q, r); }

template < class FT >
PlaneS3<FT>::PlaneS3(const PointS3<FT>& p, const DirectionS3<FT>& d)
{
  new_rep(d.dx(), d.dy(),
          d.dz(),
          -d.dx() * p.x() - d.dy() * p.y() - d.dz() * p.z());
}

template < class FT >
PlaneS3<FT>::PlaneS3(const PointS3<FT>& p, const VectorS3<FT>& v)
{ new_rep(v.x(), v.y(), v.z(), -v.x() * p.x() - v.y() * p.y() - v.z() * p.z()); }

template < class FT >
PlaneS3<FT>::PlaneS3(const FT& a, const FT& b, const FT& c, const FT& d)
{ new_rep(a, b, c, d); }

template < class FT >
PlaneS3<FT>::PlaneS3(const LineS3<FT>& l, const PointS3<FT>& p)
{ new_rep(l.point(), l.point()+l.direction().vector(), p); }

template < class FT >
PlaneS3<FT>::PlaneS3(const SegmentS3<FT>& s, const PointS3<FT>& p)
{ new_rep(s.start(), s.end(), p); }

template < class FT >
PlaneS3<FT>::PlaneS3(RayS3<FT>& r, const PointS3<FT>& p)
{ new_rep(r.start(), r.second_point(), p); }


template < class FT >
bool PlaneS3<FT>::operator==(const PlaneS3<FT>& p) const
{
  return has_on_boundary(p.point()) &&
         (orthogonal_direction() == p.orthogonal_direction());

}

template < class FT >
bool PlaneS3<FT>::operator!=(const PlaneS3<FT>& p) const
{
  return !(*this == p);
}

template < class FT >
const FT&
PlaneS3<FT>::a() const
{ return e0; }

template < class FT >
const FT&
PlaneS3<FT>::b() const
{ return e1; }

template < class FT >
const FT&
PlaneS3<FT>::c() const
{ return e2; }

template < class FT >
const FT&
PlaneS3<FT>::d() const
{ return e3; }

template < class FT >
PointS3<FT>  PlaneS3<FT>::point() const
{
  if (a() != FT(0)) // not parallel to x-axis
    return PointS3<FT>(-d()/a(), FT(0), FT(0));
  if (b() != FT(0)) // not parallel to y-axis
    return PointS3<FT>(FT(0), -d()/b(), FT(0));
  // parallel to xy-plane => intersects z-axis
  return PointS3<FT>(FT(0), FT(0), -d()/c());
}

template < class FT >
PointS3<FT>
PlaneS3<FT>::projection(const PointS3<FT>& p) const
{
  return CGAL::projection(p, *this);
}

template < class FT >
VectorS3<FT> PlaneS3<FT>::orthogonal_vector() const
{
  return VectorS3<FT>(a(), b(), c());
}


template < class FT >
DirectionS3<FT> PlaneS3<FT>::orthogonal_direction() const
{
  return DirectionS3<FT>(a(), b(), c());
}

template < class FT >
VectorS3<FT> PlaneS3<FT>::base1() const
{
  if( a() == FT(0) )  // parallel to x-axis
      return VectorS3<FT>(FT(1), FT(0), FT(0));

  if( b() == FT(0) )  // parallel to y-axis
      return VectorS3<FT>(FT(0), FT(1), FT(0));

  if (c() == FT(0) )  // parallel to z-axis
      return VectorS3<FT>(FT(0), FT(0), FT(1));

  return VectorS3<FT>(-b(), a(), FT(0));
}


template < class FT >
VectorS3<FT> PlaneS3<FT>::base2() const
{
  if ( a() == FT(0) ) // parallel to x-axis  x-axis already returned in base1
    {
      if (b() == FT(0) )  // parallel to y-axis
          return VectorS3<FT>(FT(0), FT(1), FT(0));

      if (c() == FT(0) ) // parallel to z-axis
          return VectorS3<FT>(FT(0), FT(0), FT(1));

      return VectorS3<FT>(FT(0), -b(), c());
    }
  if (b() == FT(0) )
      return VectorS3<FT>(c(), FT(0), -a());

  if (c() == FT(0) )
      return VectorS3<FT>(-b(), a(), FT(0));

  return VectorS3<FT>(FT(0), -c(), b());
}
template < class FT >
PointS3<FT> PlaneS3<FT>::to_plane_basis(const PointS3<FT>& p) const
{
  const VectorS3<FT>& v0 = base1();
  const VectorS3<FT>& v1 = base2();
  VectorS3<FT> v2 = orthogonal_vector();
  VectorS3<FT> v3 = p - point();
  FT alpha, beta, gamma;

  solve(v0.x(), v0.y(), v0.z(),
        v1.x(), v1.y(), v1.z(),
        v2.x(), v2.y(), v2.z(),
        v3.x(), v3.y(), v3.z(),
        alpha, beta, gamma);

  return PointS3<FT>(alpha, beta, gamma);
}

template < class FT >
PointS2<FT> PlaneS3<FT>::to_2d(const PointS3<FT>& p) const
{
  const VectorS3<FT>& v0 = base1();
  const VectorS3<FT>& v1 = base2();
  VectorS3<FT> v2 = orthogonal_vector();
  VectorS3<FT> v3 = p - point();
  FT alpha, beta, gamma;

  solve(v0.x(), v0.y(), v0.z(),
        v1.x(), v1.y(), v1.z(),
        v2.x(), v2.y(), v2.z(),
        v3.x(), v3.y(), v3.z(),
        alpha, beta, gamma);

  return PointS2<FT>(alpha, beta);
}


template < class FT >
PointS3<FT> PlaneS3<FT>::to_3d(const PointS2<FT>& p) const
{
  VectorS3<FT> e1 = base1(),
               e2 = base2();
  return point() + p.x() * e1 + p.y() * e2;
}

template < class FT >
LineS3<FT> PlaneS3<FT>::perpendicular_line(const PointS3<FT>& p) const
{ return LineS3<FT>(p, orthogonal_direction()); }


template < class FT >
PlaneS3<FT> PlaneS3<FT>::opposite() const
{ return PlaneS3<FT>(-a(),-b(),-c(),-d()); }


template < class FT >
PlaneS3<FT> PlaneS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{
  DirectionS3<FT> dir = t.transpose().inverse().transform(orthogonal_direction());
  if (!t.is_even())  dir = -dir;
  return PlaneS3<FT>( t.transform(point()), dir);
  
/*
  return PlaneS3<FT>( t.transform(point()), (t.is_even())
           ?   t.transpose().inverse().transform(orthogonal_direction())
           : - t.transpose().inverse().transform(orthogonal_direction()) );
*/
}


template < class FT >
Oriented_side PlaneS3<FT>::oriented_side(const PointS3<FT>& p) const
{ return Oriented_side(CGAL_NTS sign(a()*p.x() + b()*p.y() + c()*p.z() +d())); }

template < class FT >
bool PlaneS3<FT>::has_on_boundary(const  PointS3<FT>& p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) == FT(0);
}

template < class FT >
bool PlaneS3<FT>::has_on_boundary(const  LineS3<FT>& l) const
{
  return has_on_boundary(l.point())
         &&  has_on_boundary(l.point() + l.direction().vector());
}

template < class FT >
bool PlaneS3<FT>::has_on_positive_side(const  PointS3<FT>& p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) > FT(0);
}

template < class FT >
bool PlaneS3<FT>::has_on_negative_side(const  PointS3<FT>& p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) < FT(0);
}


template < class FT >
bool PlaneS3<FT>::is_degenerate() const
{
  return (a() == FT(0)) && (b() == FT(0)) && (c() == FT(0));
}


#ifndef CGAL_NO_OSTREAM_INSERT_PLANES3
template < class FT >
std::ostream& operator<<(std::ostream& os, const PlaneS3<FT>& p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.a() << ' ' << p.b() <<  ' ' << p.c() << ' ' << p.d();
    case IO::BINARY :
        write(os, p.a());
        write(os, p.b());
        write(os, p.c());
        write(os, p.d());
        return os;
        default:
            os << "PlaneS3(" << p.a() <<  ", " << p.b() <<   ", ";
            os << p.c() << ", " << p.d() <<")";
            return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANES3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANES3
template < class FT >
std::istream& operator>>(std::istream& is, PlaneS3<FT>& p)
{
    FT a, b, c, d;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c >> d;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        read(is, d);
        break;
    default:
        cerr << "" << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;
    }
    p = PlaneS3<FT>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANES3


CGAL_END_NAMESPACE

#endif  // CGAL_PLANES3_H
