#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_VECTOR_2_C
#define CGAL_CARTESIAN_VECTOR_2_C

#ifndef CGAL_CARTESIAN_DIRECTION_2_H
#include <CGAL/Cartesian/Direction_2.h>
#endif // CGAL_CARTESIAN_DIRECTION_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple<typename VectorC2<R CGAL_CTAG>::FT>*
VectorC2<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<FT>*)PTR;
}

template < class R >
VectorC2<R CGAL_CTAG>::VectorC2()
{
  PTR = new _Twotuple<typename VectorC2<R CGAL_CTAG>::FT>();
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const VectorC2<R CGAL_CTAG>  &v) :
  Handle((Handle&)v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const Null_vector &)
{
  PTR = new _Twotuple<FT>(FT(0), FT(0));
}


template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const typename VectorC2<R CGAL_CTAG>::FT &hx,
                                const typename VectorC2<R CGAL_CTAG>::FT &hy,
				const typename VectorC2<R CGAL_CTAG>::FT &hw)
{
  if ( hw == FT(1))
    PTR = new _Twotuple<FT>(hx, hy);
  else
    PTR = new _Twotuple<FT>(hx/hw, hy/hw);
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const typename VectorC2<R CGAL_CTAG>::FT &x,
                                const typename VectorC2<R CGAL_CTAG>::FT &y)
{
  PTR = new _Twotuple<FT>(x, y);
}

template < class R >
inline
VectorC2<R CGAL_CTAG>::~VectorC2()
{}

template < class R >
inline
VectorC2<R CGAL_CTAG> &
VectorC2<R CGAL_CTAG>::operator=(const VectorC2<R CGAL_CTAG> &v)
{
  Handle::operator=(v);
  return *this;
}
template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const typename VectorC2<R CGAL_CTAG>::Point_2 &p)
  : Handle((Handle&)p)
{
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC2<R CGAL_CTAG>::VectorC2(const typename VectorC2<R CGAL_CTAG>::Direction_2 &d)
  : Handle((Handle&)d)
{
}
template < class R >
CGAL_KERNEL_INLINE
bool
VectorC2<R CGAL_CTAG>::operator==(const VectorC2<R CGAL_CTAG> &v) const
{
  return (x() == v.x()) && (y() == v.y());
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator!=(const VectorC2<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator==(const Null_vector &) const
{
  return (x() == FT(0)) && (y() == FT(0));
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}

template < class R >
inline
int
VectorC2<R CGAL_CTAG>::id() const
{
  return (int) PTR ;
}
template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::x()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::y()  const
{
  return ptr()->e1;
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT 
VectorC2<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
int
VectorC2<R CGAL_CTAG>::dimension() const
{
  return 2;
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::hx()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::hy()  const
{
  return ptr()->e1;
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::hw()  const
{
  return FT(1);
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT 
VectorC2<R CGAL_CTAG>::homogeneous(int i) const
{
  return (i == 2) ? FT(1) : cartesian(i);
}
template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator+(const VectorC2<R CGAL_CTAG> &w) const
{
  return VectorC2<R CGAL_CTAG>(x() + w.x(), y() + w.y()) ;
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator-(const VectorC2<R CGAL_CTAG> &w) const
{
  return VectorC2<R CGAL_CTAG>(x() - w.x(), y() - w.y()) ;
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator-() const
{
  return VectorC2<R CGAL_CTAG>(-x(), -y()) ;
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::operator*(const VectorC2<R CGAL_CTAG> &w) const
{
  return x() * w.x() + y() * w.y() ;
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator/(const typename VectorC2<R CGAL_CTAG>::FT &c) const
{
  return VectorC2<R CGAL_CTAG>( x()/c, y()/c) ;
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::Direction_2
VectorC2<R CGAL_CTAG>::direction() const
{
  return Direction_2(*this) ;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::perpendicular(const Orientation &o) const
{
  CGAL_kernel_precondition( o != COLLINEAR );
  if (o == COUNTERCLOCKWISE)
    return VectorC2<R CGAL_CTAG>(-y(), x());
  else
    return VectorC2<R CGAL_CTAG>(y(), -x());
}

template < class R >
inline
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::transform(const typename VectorC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return t.transform(*this);
}



#ifndef CGAL_NO_OSTREAM_INSERT_VECTORC2
template < class R >
std::ostream &operator<<(std::ostream &os, const VectorC2<R CGAL_CTAG> &v)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "VectorC2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORC2
template < class R >
std::istream &operator>>(std::istream &is, VectorC2<R CGAL_CTAG> &p)
{
    typename VectorC2<R CGAL_CTAG>::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        cerr << "" << std::endl;
        cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = VectorC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_VECTOR_2_C
