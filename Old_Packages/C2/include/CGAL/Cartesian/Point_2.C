#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef _MSC_VER
#define typename
#endif

#ifndef CGAL_CARTESIAN_POINT_2_C
#define CGAL_CARTESIAN_POINT_2_C

#ifndef CGAL_ORIGIN_H
#include <CGAL/Origin.h>
#endif // CGAL_ORIGIN_H
#ifndef CGAL_BBOX_2_H
#include <CGAL/Bbox_2.h>
#endif // CGAL_BBOX_2_H
#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif // CGAL_NUMBER_UTILS_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple<typename PointC2<R CGAL_CTAG>::FT>*
PointC2<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<FT>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2()
{
  PTR = new _Twotuple<FT>;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2(const Origin &)
{
  PTR = new _Twotuple<FT>(FT(0), FT(0));
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2(const PointC2<R CGAL_CTAG> &p)
  : Handle((Handle&)p)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2(const typename PointC2<R CGAL_CTAG>::Vector_2 &v)
  : Handle((Handle&)v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2(const typename PointC2<R CGAL_CTAG>::FT &hx,
                    const typename PointC2<R CGAL_CTAG>::FT &hy,
                    const typename PointC2<R CGAL_CTAG>::FT &hw)
{
  if( hw != FT(1)) {
    PTR = new _Twotuple<FT>(hx/hw, hy/hw);
  } else {
    PTR = new _Twotuple<FT>(hx, hy);
  }
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointC2<R CGAL_CTAG>::PointC2(const typename PointC2<R CGAL_CTAG>::FT &x,
                              const typename PointC2<R CGAL_CTAG>::FT &y)
{
  PTR = new _Twotuple<FT>(x, y);
}

template < class R >
inline
PointC2<R CGAL_CTAG>::~PointC2()
{}

template < class R >
inline
PointC2<R CGAL_CTAG> &
PointC2<R CGAL_CTAG>::operator=(const PointC2<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
inline
bool PointC2<R CGAL_CTAG>::operator==(const PointC2<R CGAL_CTAG>& p) const
{
  return ((x() == p.x()) && (y() == p.y())) ;
}

template < class R >
inline
bool PointC2<R CGAL_CTAG>::operator!=(const PointC2<R CGAL_CTAG>& p) const
{
  return !(*this == p);
}

template < class R >
inline
int PointC2<R CGAL_CTAG>::id() const { return (int)PTR; }

template < class R >
inline typename PointC2<R CGAL_CTAG>::FT PointC2<R CGAL_CTAG>::x()  const
{
  return ptr()->e0;
}

template < class R >
inline typename PointC2<R CGAL_CTAG>::FT PointC2<R CGAL_CTAG>::y()  const
{
  return  ptr()->e1 ;
}

template < class R >
CGAL_KERNEL_INLINE
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
inline
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
int
PointC2<R CGAL_CTAG>::dimension() const
{
  return 2;
}

template < class R >
inline
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::hx()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::hy()  const
{
  return ptr()->e1;
}

template < class R >
inline
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::hw()  const
{
  return FT(1);
}

template < class R >
inline
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2)
    return cartesian(i);
  return FT(1);
}

template < class R >
CGAL_KERNEL_INLINE
PointC2<R CGAL_CTAG>
PointC2<R CGAL_CTAG>::
transform( const typename PointC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return t.transform(*this);
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2 PointC2<R CGAL_CTAG>::bbox() const
{
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  return Bbox_2(bx,by, bx,by);
}

#ifndef CGAL_NO_OSTREAM_INSERT_POINTC2
template < class R >
std::ostream &operator<<(std::ostream &os, const PointC2<R CGAL_CTAG> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        return os;
    default:
        return os << "PointC2(" << p.x() << ", " << p.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTC2
template < class R >
std::istream &operator>>(std::istream &is, PointC2<R CGAL_CTAG> &p)
{
    typename PointC2<R CGAL_CTAG>::FT x, y;
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
    p = PointC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTC2

CGAL_END_NAMESPACE

#ifdef _MSC_VER
#undef typename
#endif

#endif // CGAL_CARTESIAN_POINT_2_C
