#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_DIRECTION_2_C
#define CGAL_CARTESIAN_DIRECTION_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple<typename R::FT>*
DirectionC2<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<FT>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::DirectionC2()
{
  PTR = new _Twotuple<FT>();
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::
DirectionC2(const DirectionC2<R CGAL_CTAG> &d)
  : Handle((Handle&)d)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::
DirectionC2(const DirectionC2<R CGAL_CTAG>::Vector_2 &v) :
  Handle((Handle&)v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::DirectionC2(const typename R::FT &x,
                                      const typename R::FT &y)
{
  PTR = new _Twotuple<FT>(x, y);
}

template < class R >
inline DirectionC2<R CGAL_CTAG>:: ~DirectionC2()
{}

template < class R >
CGAL_KERNEL_INLINE
DirectionC2<R CGAL_CTAG> &
DirectionC2<R CGAL_CTAG>::operator=(const DirectionC2<R CGAL_CTAG> &d)
{
  Handle::operator=(d);
  return *this;
}
template < class R >
bool
DirectionC2<R CGAL_CTAG>::operator==(const DirectionC2<R CGAL_CTAG> &d) const
{
// Use a C2 predicate for that ?
  return (CGAL::sign(dx()) == CGAL::sign(d.dx()))
      && (CGAL::sign(dy()) == CGAL::sign(d.dy()))
      && (dy()*d.dx() == d.dy()*dx());
}

template < class R >
inline
bool
DirectionC2<R CGAL_CTAG>::operator!=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return !( *this == d );
}

template < class R >
inline
int
DirectionC2<R CGAL_CTAG>::id() const
{
  return (int)PTR;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator<(const DirectionC2<R CGAL_CTAG> &d) const
{
  int quadrant_this = (dx() >= FT(0)) ? ((dy() >= FT(0))?1:4)
                                      : ((dy() >= FT(0))?2:3);
  int quadrant_d    = (d.dx() >= FT(0)) ? ((d.dy() >= FT(0))?1:4)
                                        : ((d.dy() >= FT(0))?2:3);

  if(quadrant_this < quadrant_d)
    return true;
  else if (quadrant_this > quadrant_d)
    return false;
  else
    return dy() * d.dx() < d.dy() * dx();
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator>(const DirectionC2<R CGAL_CTAG> &d) const
{
  return d < *this ;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator>=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return (d < *this) || (d == *this) ;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator<=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return (*this < d) || (d == *this) ;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::counterclockwise_in_between
   (const DirectionC2<R CGAL_CTAG> &d1,
    const DirectionC2<R CGAL_CTAG> &d2) const
{
  return (d2 > *this) && (*this > d1) ;
}

template < class R >
inline
DirectionC2<R CGAL_CTAG>::Vector_2
DirectionC2<R CGAL_CTAG>::vector() const
{
  return Vector_2(*this);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
DirectionC2<R CGAL_CTAG>
DirectionC2<R CGAL_CTAG>::perpendicular(const Orientation &o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
    return DirectionC2<R CGAL_CTAG>(-dy(), dx());
  else
    return DirectionC2<R CGAL_CTAG>(dy(), -dx());
}

template < class R >
CGAL_KERNEL_INLINE
DirectionC2<R CGAL_CTAG>
DirectionC2<R CGAL_CTAG>::
transform(const DirectionC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return t.transform(*this);
}

template < class R >
inline
DirectionC2<R CGAL_CTAG>
DirectionC2<R CGAL_CTAG>::operator-() const
{
  return DirectionC2<R CGAL_CTAG>(-dx(), -dy());
}



template < class R >
CGAL_KERNEL_INLINE
typename R::FT
DirectionC2<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  return (i==0) ? dx() : dy();
}


template < class R >
inline
typename R::FT
DirectionC2<R CGAL_CTAG>::dx() const
{
  return ptr()->e0;
}

template < class R >
inline
typename R::FT
DirectionC2<R CGAL_CTAG>::dy() const
{
  return ptr()->e1;
}


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC2
template < class R >
std::ostream
&operator<<(std::ostream &os, const DirectionC2<R CGAL_CTAG> &d)
{
    typename R::Vector_2 v = d.vector();
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
std::istream
&operator>>(std::istream &is, DirectionC2<R CGAL_CTAG> &p)
{
    typename R::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        cerr << std::endl << "Stream must be in ascii or binary mode"
	     << std::endl;
        break;
    }
    p = DirectionC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DIRECTION_2_C
