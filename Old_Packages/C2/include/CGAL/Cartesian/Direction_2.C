// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_DIRECTION_2_C
#define CGAL_CARTESIAN_DIRECTION_2_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE


template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::DirectionC2()
{
  new ( static_cast< void*>(ptr)) Twotuple<FT>(); 
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::
DirectionC2(const DirectionC2<R CGAL_CTAG> &d)
  : Handle_for<Twotuple<typename R::FT> >(d)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::
DirectionC2(const typename DirectionC2<R CGAL_CTAG>::Vector_2 &v)
 : Handle_for<Twotuple<typename R::FT> >(v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC2<R CGAL_CTAG>::
DirectionC2(const typename DirectionC2<R CGAL_CTAG>::FT &x,
            const typename DirectionC2<R CGAL_CTAG>::FT &y)
{
  new ( static_cast< void*>(ptr)) Twotuple<FT>(x, y);
}

template < class R >
inline
DirectionC2<R CGAL_CTAG>::~DirectionC2()
{}


template < class R >
inline
bool
DirectionC2<R CGAL_CTAG>::operator==(const DirectionC2<R CGAL_CTAG> &d) const
{
  if ( ptr == d.ptr ) return true;
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
  return compare_angle_with_x_axis(*this,d) == SMALLER;
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
  return compare_angle_with_x_axis(*this,d) != SMALLER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R CGAL_CTAG>::operator<=(const DirectionC2<R CGAL_CTAG> &d) const
{
  return compare_angle_with_x_axis(*this,d) != LARGER;
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
  if ( d1 < *this)  {
    return ( *this < d2 )||( d2 <= d1 );
  }
  else    {
    return ( *this < d2 )&&( d2 <= d1 );
  }
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
transform(const typename DirectionC2<R CGAL_CTAG>::Aff_transformation_2 &t)
    const
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
typename DirectionC2<R CGAL_CTAG>::FT
DirectionC2<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  return (i==0) ? dx() : dy();
}

template < class R >
inline
typename DirectionC2<R CGAL_CTAG>::FT
DirectionC2<R CGAL_CTAG>::dx() const
{
  return ptr->e0;
}

template < class R >
inline
typename DirectionC2<R CGAL_CTAG>::FT
DirectionC2<R CGAL_CTAG>::dy() const
{
  return ptr->e1;
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
    p = DirectionC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_DIRECTION_2_C
