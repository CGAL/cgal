#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_PLANE_3_C
#define CGAL_CARTESIAN_PLANE_3_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Fourtuple<typename R::FT>*
PlaneC3<R CGAL_CTAG>::ptr() const
{
    return (_Fourtuple<FT>*)PTR;
}

template < class R >
inline
void
PlaneC3<R CGAL_CTAG>::
new_rep(const typename R::FT &a, const typename R::FT &b,
        const typename R::FT &c, const typename R::FT &d)
{
  PTR = new _Fourtuple<FT>(a, b, c, d);
}

template < class R >
inline
void
PlaneC3<R CGAL_CTAG>::new_rep(const PlaneC3<R CGAL_CTAG>::Point_3 &p,
                              const PlaneC3<R CGAL_CTAG>::Point_3 &q,
                              const PlaneC3<R CGAL_CTAG>::Point_3 &r)
{
  FT rpx = p.x()-r.x();
  FT rpy = p.y()-r.y();
  FT rpz = p.z()-r.z();
  FT rqx = q.x()-r.x();
  FT rqy = q.y()-r.y();
  FT rqz = q.z()-r.z();
  // Cross product rp * rq.
  FT A = rpy*rqz - rqy*rpz;
  FT B = rpz*rqx - rqz*rpx;
  FT C = rpx*rqy - rqx*rpy;
  FT D = - A*r.x() - B*r.y() - C*r.z();
  PTR = new _Fourtuple<FT>(A, B, C, D);
}


CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_LINE_3_H
#include <CGALCartesian//Line_3.h>
#endif // CGAL_CARTESIAN_LINE_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3()
{
  PTR = new _Fourtuple<FT>();
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG> &p)
  : Handle(p)
{}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const PlaneC3<R CGAL_CTAG>::Point_3 &q,
        const PlaneC3<R CGAL_CTAG>::Point_3 &r)
{
  new_rep(p, q, r);
}

template < class R >
CGAL_KERNEL_INLINE
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const PlaneC3<R CGAL_CTAG>::Direction_3 &d)
{
  new_rep(d.dx(), d.dy(),
          d.dz(),
          -d.dx() * p.x() - d.dy() * p.y() - d.dz() * p.z());
}

template < class R >
CGAL_KERNEL_INLINE
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const PlaneC3<R CGAL_CTAG>::Vector_3 &v)
{
  new_rep(v.x(), v.y(), v.z(), -v.x() * p.x() - v.y() * p.y() - v.z() * p.z());
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename R::FT &a,
        const typename R::FT &b,
        const typename R::FT &c,
        const typename R::FT &d)
{
  new_rep(a, b, c, d);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Line_3 &l,
        const PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(l.point(), l.point()+l.direction().vector(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Segment_3 &s,
        const PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(s.start(), s.end(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG>::Ray_3 &r,
        const PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(r.start(), r.second_point(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::~PlaneC3()
{}

template < class R >
inline
PlaneC3<R CGAL_CTAG> &PlaneC3<R CGAL_CTAG>::operator=(const PlaneC3<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool PlaneC3<R CGAL_CTAG>::operator==(const PlaneC3<R CGAL_CTAG> &p) const
{
  return has_on_boundary(p.point()) &&
         (orthogonal_direction() == p.orthogonal_direction());

}

template < class R >
inline
bool PlaneC3<R CGAL_CTAG>::operator!=(const PlaneC3<R CGAL_CTAG> &p) const
{
  return !(*this == p);
}

template < class R >
inline
long PlaneC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename R::FT PlaneC3<R CGAL_CTAG>::a() const
{
  return ptr()->e0;
}

template < class R >
inline
typename R::FT PlaneC3<R CGAL_CTAG>::b() const
{
  return ptr()->e1;
}

template < class R >
inline
typename R::FT PlaneC3<R CGAL_CTAG>::c() const
{
  return ptr()->e2;
}

template < class R >
inline
typename R::FT PlaneC3<R CGAL_CTAG>::d() const
{
  return ptr()->e3;
}

template < class R >
PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::point() const
{
  if (a() != FT(0)) // not parallel to x-axis
    return PlaneC3<R CGAL_CTAG>::Point_3(-d()/a(), FT(0), FT(0));
  if (b() != FT(0)) // not parallel to y-axis
    return PlaneC3<R CGAL_CTAG>::Point_3(FT(0), -d()/b(), FT(0));
  // parallel to xy-plane => intersects z-axis
  return PlaneC3<R CGAL_CTAG>::Point_3(FT(0), FT(0), -d()/c());
}

template < class R >
PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
projection(const PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return CGAL::projection(p, *this);
}

template < class R >
PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::orthogonal_vector() const
{
  return PlaneC3<R CGAL_CTAG>::Vector_3(a(), b(), c());
}

template < class R >
PlaneC3<R CGAL_CTAG>::Direction_3
PlaneC3<R CGAL_CTAG>::orthogonal_direction() const
{
  return PlaneC3<R CGAL_CTAG>::Direction_3(a(), b(), c());
}

template < class R >
PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::base1() const
{
  if( a() == FT(0) )  // parallel to x-axis
      return Vector_3(FT(1), FT(0), FT(0));

  if( b() == FT(0) )  // parallel to y-axis
      return Vector_3(FT(0), FT(1), FT(0));

  if (c() == FT(0) )  // parallel to z-axis
      return Vector_3(FT(0), FT(0), FT(1));

  return PlaneC3<R CGAL_CTAG>::Vector_3(-b(), a(), FT(0));
}

template < class R >
PlaneC3<R CGAL_CTAG>::Vector_3 PlaneC3<R CGAL_CTAG>::base2() const
{
  if ( a() == FT(0) ) // parallel to x-axis  x-axis already returned in base1
    {
      if (b() == FT(0) )  // parallel to y-axis
          return Vector_3(FT(0), FT(1), FT(0));

      if (c() == FT(0) ) // parallel to z-axis
          return Vector_3(FT(0), FT(0), FT(1));

      return Vector_3(FT(0), -b(), c());
    }
  if (b() == FT(0) )
      return Vector_3(c(), FT(0), -a());

  if (c() == FT(0) )
      return Vector_3(-b(), a(), FT(0));

  return Vector_3(FT(0), -c(), b());
}

template < class R >
PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_plane_basis(const PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  const Vector_3 &e1 = base1();
  const Vector_3 &e2 = base2();
  FT alpha, beta, gamma;

  solve(e1, e2, orthogonal_vector(), p - point(), alpha, beta, gamma);

  return Point_3(alpha, beta, gamma);
}

template < class R >
PlaneC3<R CGAL_CTAG>::Point_2
PlaneC3<R CGAL_CTAG>::
to_2d(const PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  const PlaneC3<R CGAL_CTAG>::Vector_3 &e1 = base1();
  const PlaneC3<R CGAL_CTAG>::Vector_3 &e2 = base2();
  FT alpha, beta, gamma;

  solve(e1, e2, orthogonal_vector(), p - point(), alpha, beta, gamma);

  return Point_2(alpha, beta);
}

template < class R >
PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_3d(const PlaneC3<R CGAL_CTAG>::Point_2 &p) const
{
  PlaneC3<R CGAL_CTAG>::Vector_3 e1 = base1(),
               e2 = base2();
  return point() + p.x() * e1 + p.y() * e2;
}

template < class R >
PlaneC3<R CGAL_CTAG>::Line_3
PlaneC3<R CGAL_CTAG>::
perpendicular_line(const PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return PlaneC3<R CGAL_CTAG>::Line_3(p, orthogonal_direction());
}

template < class R >
PlaneC3<R CGAL_CTAG>
PlaneC3<R CGAL_CTAG>::opposite() const
{
  return PlaneC3<R CGAL_CTAG>(-a(),-b(),-c(),-d());
}

template < class R >
PlaneC3<R CGAL_CTAG>
PlaneC3<R CGAL_CTAG>::
transform(const PlaneC3<R CGAL_CTAG>::Aff_transformation_3& t) const
{
  return PlaneC3<R CGAL_CTAG>( t.transform(point()), (t.is_even())
           ?   t.transpose().inverse().transform(orthogonal_direction())
           : - t.transpose().inverse().transform(orthogonal_direction()) );
}

template < class R >
Oriented_side
PlaneC3<R CGAL_CTAG>::
oriented_side(const PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Oriented_side(CGAL::sign(a()*p.x() + b()*p.y() + c()*p.z() +d()));
}

template < class R >
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) == FT(0);
}

template < class R >
bool
PlaneC3<R CGAL_CTAG>::
has_on(const  PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return has_on_boundary(p);
}

template < class R >
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  PlaneC3<R CGAL_CTAG>::Line_3 &l) const
{
  return has_on_boundary(l.point())
         &&  has_on_boundary(l.point() + l.direction().vector());
}

template < class R >
bool
PlaneC3<R CGAL_CTAG>::
has_on_positive_side(const  PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) > FT(0);
}

template < class R >
bool
PlaneC3<R CGAL_CTAG>::
has_on_negative_side(const  PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return (a()*p.x() + b()*p.y() + c()*p.z() +d()) < FT(0);
}

template < class R >
bool PlaneC3<R CGAL_CTAG>::is_degenerate() const
{
  return (a() == FT(0)) && (b() == FT(0)) && (c() == FT(0));
}


#ifndef CGAL_NO_OSTREAM_INSERT_PLANEC3
template < class R >
std::ostream &operator<<(std::ostream &os, const PlaneC3<R CGAL_CTAG> &p)
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
            os << "PlaneC3(" << p.a() <<  ", " << p.b() <<   ", ";
            os << p.c() << ", " << p.d() <<")";
            return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANEC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANEC3
template < class R >
std::istream &operator>>(std::istream &is, PlaneC3<R CGAL_CTAG> &p)
{
    typename R::FT a, b, c, d;
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
        cerr << "" << std::endl;
        cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = PlaneC3<R CGAL_CTAG>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANEC3


CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PLANE_3_C
