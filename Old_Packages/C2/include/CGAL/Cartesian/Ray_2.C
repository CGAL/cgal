#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_RAY_2_C
#define CGAL_CARTESIAN_RAY_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename RayC2<R CGAL_CTAG>::Point_2 > *
RayC2<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<Point_2>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
RayC2<R CGAL_CTAG>::RayC2()
{
  PTR = new _Twotuple<Point_2>;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
RayC2<R CGAL_CTAG>::RayC2(const RayC2<R CGAL_CTAG>  &r)
  : Handle((Handle&)r)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
RayC2<R CGAL_CTAG>::RayC2(const typename RayC2<R CGAL_CTAG>::Point_2 &sp,
                          const typename RayC2<R CGAL_CTAG>::Point_2 &secondp)
{
  PTR = new _Twotuple<Point_2>(sp, secondp);
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
RayC2<R CGAL_CTAG>::RayC2(const typename RayC2<R CGAL_CTAG>::Point_2 &sp,
                          const typename RayC2<R CGAL_CTAG>::Direction_2 &d)
{
  PTR = new _Twotuple<Point_2>(sp, sp + d.vector());
}

template < class R >
inline
RayC2<R CGAL_CTAG>::~RayC2()
{}

template < class R >
inline
RayC2<R CGAL_CTAG> &RayC2<R CGAL_CTAG>::operator=(const RayC2<R CGAL_CTAG> &r)
{
  Handle::operator=(r);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R CGAL_CTAG>::operator==(const RayC2<R CGAL_CTAG> &r) const
{
  return ((source() == r.source()) && (direction() == r.direction()) );
}

template < class R >
bool RayC2<R CGAL_CTAG>::operator!=(const RayC2<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
int RayC2<R CGAL_CTAG>::id() const
{
  return (int) PTR ;
}

template < class R >
inline
typename RayC2<R CGAL_CTAG>::Point_2
RayC2<R CGAL_CTAG>::start() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayC2<R CGAL_CTAG>::Point_2
RayC2<R CGAL_CTAG>::source() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayC2<R CGAL_CTAG>::Point_2
RayC2<R CGAL_CTAG>::second_point() const
{
  return ptr()->e1;
}

template < class R >
CGAL_KERNEL_INLINE
typename RayC2<R CGAL_CTAG>::Point_2
RayC2<R CGAL_CTAG>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0) return ptr()->e0;
  if (i == 1) return ptr()->e1;
  return source() + FT(i) * (second_point() - source());
}

template < class R >
inline
typename RayC2<R CGAL_CTAG>::Direction_2
RayC2<R CGAL_CTAG>::direction() const
{
  return Direction_2( second_point() - source() );
}

template < class R >
inline
typename RayC2<R CGAL_CTAG>::Line_2
RayC2<R CGAL_CTAG>::supporting_line() const
{
  return Line_2(*this);
}

template < class R >
inline
RayC2<R CGAL_CTAG>
RayC2<R CGAL_CTAG>::opposite() const
{
  return RayC2<R CGAL_CTAG>( source(), - direction() );
}


template < class R >
CGAL_KERNEL_INLINE
RayC2<R CGAL_CTAG>
RayC2<R CGAL_CTAG>::
transform(const typename RayC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return RayC2<R CGAL_CTAG>(t.transform(source()), t.transform(second_point()));
}


template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R CGAL_CTAG>::is_horizontal() const
{
  return (source().y() ==  second_point().y());
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R CGAL_CTAG>::is_vertical() const
{
  return  (source().x() == second_point().x());
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R CGAL_CTAG>::is_degenerate() const
{
  return (source() == second_point());
}

template < class R >
CGAL_KERNEL_INLINE
bool
RayC2<R CGAL_CTAG>::has_on(const typename RayC2<R CGAL_CTAG>::Point_2 &p) const
{
  return ( p == source()
           || ( collinear(source(), p, second_point())
           && ( typename R::Direction_2(p - source()) == direction() )));

}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
RayC2<R CGAL_CTAG>::collinear_has_on(const typename RayC2<R CGAL_CTAG>::Point_2 &p) const
{
    switch(compare_x(source(), second_point())){
    case SMALLER:
        return compare_x(source(), p) != LARGER;
    case LARGER:
        return compare_x(p, source()) != LARGER;
    default:
        switch(compare_y(source(), second_point())){
        case SMALLER:
            return compare_y(source(), p) != LARGER;
        case LARGER:
            return compare_y(p, source()) != LARGER;
        default:
            return true; // p == source()
        }
    }
}

#ifndef CGAL_NO_OSTREAM_INSERT_RAYC2
template < class R >
std::ostream &operator<<(std::ostream &os, const RayC2<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.source() << ' ' << r.direction();
    case IO::BINARY :
        return os << r.source() << r.direction();
    default:
        return os << "RayC2(" << r.source() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC2
template < class R >
std::istream &operator>>(std::istream &is, RayC2<R CGAL_CTAG> &r)
{
    typename RayC2<R CGAL_CTAG>::Point_2 p;
    typename RayC2<R CGAL_CTAG>::Direction_2 d;

    is >> p >> d;

    r = RayC2<R CGAL_CTAG>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC2

CGAL_END_NAMESPACE

#endif
