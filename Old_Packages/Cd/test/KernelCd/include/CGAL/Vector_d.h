#ifndef CGAL_VECTOR_D_H
#define CGAL_VECTOR_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Point_d.h>
#include <CGAL/Direction_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Vector_d : public _R::Vector_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef Vector_d<R>                       Self;
  typedef typename R::Vector_d_base         RVector_d;
  typedef const FT*                         const_iterator;

  Vector_d()
    {}
  Vector_d(const Self& v) : RVector_d( (const RVector_d& )v )
    {}
  Vector_d(const RVector_d&  v) : RVector_d(v)
    {}
  Vector_d(int dim, const Null_vector& v) : RVector_d(dim,v)
    {}
  template <class InputIterator>
  Vector_d (int dim, InputIterator first, InputIterator last)
      : RVector_d (dim, first, last)
    {}
  //Vector_d(int dim, const RT &x, const RT &y)
    //{ CGAL_kernel_precondition( dim == 2);
      //RT e[2] = { x, y }; *this = RVector_d(2,e+0,e+2);
    //}
  Vector_d(int dim, const RT &x, const RT &y, const RT &z)
    { CGAL_kernel_precondition( dim == 2 || dim == 3);
      RT e[3] = { x, y, z }; *this = RVector_d(dim,e+0,e+3);
    }
  Vector_d(int dim, const RT &x, const RT &y, const RT &z, const RT &w)
    { CGAL_kernel_precondition( dim == 3 || dim == 4);
      RT e[4] = { x, y, z, w }; *this = RVector_d(dim,e+0,e+4);
    }
  Vector_d(int dim,
           const RT &x, const RT &y, const RT &z, const RT &t, const RT &w)
    { CGAL_kernel_precondition( dim == 4);
      RT e[5] = { x, y, z, t, w }; *this = RVector_d(dim,e+0,e+5);
    }
   Self& operator=(const Self& v)
    { RVector_d::operator=(v); return *this; }
  bool operator==(const Self& v) const
    { return RVector_d::operator==(v); }
  bool operator!=(const Self& v) const
    { return !(*this == v); }
  bool operator==(const Null_vector& v) const
    { return RVector_d::operator==(v); }
  bool operator!=(const Null_vector& v) const
    { return !(*this == v); }
  RT homogeneous(int i) const
    { return RVector_d::homogeneous(i); }
  FT cartesian(int i) const
    { return RVector_d::cartesian(i); }
  FT x() { return cartesian(0); }
  FT y() { return cartesian(1); }
  FT z() { return cartesian(2); }
  RT hx() { return homogeneous(0); }
  RT hy() { return homogeneous(1); }
  RT hz() { return homogeneous(2); }
  RT hw() { return homogeneous(3); }
  FT operator[](int i) const
    { return cartesian(i); }
  int dimension() const
    { return RVector_d::dimension(); }
  Self operator+(const Self& w) const
    { return (const RVector_d& )(*this) + (const RVector_d& )(w); }
  Self operator-(const Self& w) const
    { return (const RVector_d& )(*this) - (const RVector_d& )(w); }
  Self operator-() const
    { return RVector_d::operator-(); }
  FT operator*(const Self& w) const
    { return (const RVector_d& )(*this) * (const RVector_d& )(w); }
  Self operator*(const RT& c) const
    { return (const RVector_d& )(*this) * c; }
  Self operator/(const RT& c) const
    { return (const RVector_d& )(*this) / c; }
  CGAL::Direction_d<R> direction() const
    { return RVector_d::direction(); }
  // Self transform(const CGAL::Aff_transformation_d<R>& t) const
  //  { return RVector_d::transform(t); }

private:
  Vector_d(const CGAL::Point_d<R>& p) : RVector_d(p)
    {}
  Vector_d(const CGAL::Direction_d<R>& d) : RVector_d(d)
    {}
};

template < class R >
No_number_tag number_type_tag(const CGAL::Vector_d<R>& )
{
  return No_number_tag();
}

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Vector_d<R>& v)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return os << (const RVector_d& )v;
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_D

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_D
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Vector_d<R>& p)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return is >> (RVector_d& )p;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_D

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_D_H
