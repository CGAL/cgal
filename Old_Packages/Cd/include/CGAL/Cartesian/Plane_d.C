// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann

#ifndef CGAL_CARTESIAN_PLANE_D_C
#define CGAL_CARTESIAN_PLANE_D_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#include <CGAL/Cartesian/constructions_on_planes_d.h>
#include <CGAL/Cartesian/distance_computations_d.h>
#include <CGAL/Cartesian/predicates_on_planes_d.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Fourtuple<typename R::FT>*
PlaneCd<R CGAL_CTAG>::ptr() const
{
    return (_Fourtuple<FT>*)PTR;
}

template < class R >
inline
void
PlaneCd<R CGAL_CTAG>::
new_rep(const typename PlaneCd<R CGAL_CTAG>::FT &a,
        const typename PlaneCd<R CGAL_CTAG>::FT &b,
        const typename PlaneCd<R CGAL_CTAG>::FT &c,
        const typename PlaneCd<R CGAL_CTAG>::FT &d)
{
  PTR = new _Fourtuple<FT>(a, b, c, d);
}

template < class R >
inline
void
PlaneCd<R CGAL_CTAG>::new_rep(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
                              const typename PlaneCd<R CGAL_CTAG>::Point_d &q,
                              const typename PlaneCd<R CGAL_CTAG>::Point_d &r)
{
  PlaneCd<R CGAL_CTAG> h = plane_from_points(p,q,r);
  new_rep(h.a(), h.b(), h.c(), h.d());
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd()
{
  PTR = new _Fourtuple<FT>();
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const PlaneCd<R CGAL_CTAG> &p)
  : Handle(p)
{}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
        const typename PlaneCd<R CGAL_CTAG>::Point_d &q,
        const typename PlaneCd<R CGAL_CTAG>::Point_d &r)
{
  new_rep(p, q, r);
}

template < class R >
CGAL_KERNEL_INLINE
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
        const typename PlaneCd<R CGAL_CTAG>::Direction_d &d)
{
  PlaneCd<R CGAL_CTAG> h = plane_from_point_direction(p,d);
  new_rep(h.a(), h.b(), h.c(), h.d());
}

template < class R >
CGAL_KERNEL_INLINE
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
        const typename PlaneCd<R CGAL_CTAG>::Vector_d &v)
{
  FT a, b, c, d;
  plane_from_point_directionCd(p.x(),p.y(),p.z(),v.x(),v.y(),v.z(),a,b,c,d);
  new_rep(a, b, c, d);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::FT &a,
        const typename PlaneCd<R CGAL_CTAG>::FT &b,
        const typename PlaneCd<R CGAL_CTAG>::FT &c,
        const typename PlaneCd<R CGAL_CTAG>::FT &d)
{
  new_rep(a, b, c, d);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Line_d &l,
        const typename PlaneCd<R CGAL_CTAG>::Point_d &p)
{
  new_rep(l.point(), l.point()+l.direction().to_vector(), p);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Segment_d &s,
        const typename PlaneCd<R CGAL_CTAG>::Point_d &p)
{
  new_rep(s.start(), s.end(), p);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Ray_d &r,
        const typename PlaneCd<R CGAL_CTAG>::Point_d &p)
{
  new_rep(r.start(), r.second_point(), p);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::~PlaneCd()
{}

template < class R >
inline
PlaneCd<R CGAL_CTAG> &PlaneCd<R CGAL_CTAG>::
operator=(const PlaneCd<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool PlaneCd<R CGAL_CTAG>::
operator==(const PlaneCd<R CGAL_CTAG> &p) const
{
  return has_on_boundary(p.point()) &&
         (orthogonal_direction() == p.orthogonal_direction());

}

template < class R >
inline
bool PlaneCd<R CGAL_CTAG>::
operator!=(const PlaneCd<R CGAL_CTAG> &p) const
{
  return !(*this == p);
}

template < class R >
inline
long PlaneCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::FT
PlaneCd<R CGAL_CTAG>::a() const
{
  return ptr()->e0;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::FT
PlaneCd<R CGAL_CTAG>::b() const
{
  return ptr()->e1;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::FT
PlaneCd<R CGAL_CTAG>::c() const
{
  return ptr()->e2;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::FT
PlaneCd<R CGAL_CTAG>::d() const
{
  return ptr()->d;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Point_d
PlaneCd<R CGAL_CTAG>::point() const
{
  return point_on_plane(*this);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Point_d
PlaneCd<R CGAL_CTAG>::
projection(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return CGAL::projection_plane(p, *this);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Vector_d
PlaneCd<R CGAL_CTAG>::orthogonal_vector() const
{
  return Vector_d(a(), b(), c());
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Direction_d
PlaneCd<R CGAL_CTAG>::orthogonal_direction() const
{
  return Direction_d(a(), b(), c());
}

template < class R >
typename PlaneCd<R CGAL_CTAG>::Vector_d
PlaneCd<R CGAL_CTAG>::base(const int i) const
{
  if( a() == FT(0) )  // parallel to x-axis
      return Vector_d(FT(1), FT(0), FT(0));

  if( b() == FT(0) )  // parallel to y-axis
      return Vector_d(FT(0), FT(1), FT(0));

  if (c() == FT(0) )  // parallel to z-axis
      return Vector_d(FT(0), FT(0), FT(1));

  return Vector_d(-b(), a(), FT(0));
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Line_d
PlaneCd<R CGAL_CTAG>::
perpendicular_line(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return Line_d(p, orthogonal_direction());
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>
PlaneCd<R CGAL_CTAG>::opposite() const
{
  return PlaneCd<R CGAL_CTAG>(-a(),-b(),-c(),-d());
}

template < class R >
PlaneCd<R CGAL_CTAG>
PlaneCd<R CGAL_CTAG>::
transform(const typename PlaneCd<R CGAL_CTAG>::Aff_transformation_d& t) const
{
  return PlaneCd<R CGAL_CTAG>( t.transform(point()), (t.is_even())
           ?   t.transpose().inverse().transform(orthogonal_direction())
           : - t.transpose().inverse().transform(orthogonal_direction()) );
}

template < class R >
inline
Oriented_side
PlaneCd<R CGAL_CTAG>::
oriented_side(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return side_of_oriented_plane(*this,p);
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return has_on_boundary(p);
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneCd<R CGAL_CTAG>::Line_d &l) const
{
  return has_on_boundary(l.point())
     &&  has_on_boundary(l.point() + l.direction().to_vector());
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_positive_side(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_negative_side(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
is_degenerate() const
{
  return (a() == FT(0)) && (b() == FT(0)) && (c() == FT(0));
}


#ifndef CGAL_NO_OSTREAM_INSERT_PLANECD
template < class R >
std::ostream &operator<<(std::ostream &os, const PlaneCd<R CGAL_CTAG> &p)
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
            os << "PlaneCd(" << p.a() <<  ", " << p.b() <<   ", ";
            os << p.c() << ", " << p.d() <<")";
            return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANECD

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANECD
template < class R >
std::istream &operator>>(std::istream &is, PlaneCd<R CGAL_CTAG> &p)
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = PlaneCd<R CGAL_CTAG>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANECD


CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_PLANE_D_C
