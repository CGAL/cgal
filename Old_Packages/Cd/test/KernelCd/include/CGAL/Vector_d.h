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
  typedef CGAL::Vector_3<R>                 Self;
  typedef typename R::Vector_d_base         RVector_d;
  typedef const FT*                         const_iterator;

  Vector_d()
    {}
  Vector_d(const Self& v) : RVector_d( (const RVector_d& )v )
    {}
  Vector_d(const RVector_d&  v) : RVector_d(v)
    {}
  Vector_d(const Null_vector& v) : RVector_d(v)
    {}
  template <class InputIterator>
  Vector_d (int dim, InputIterator first, InputIterator last)
      : RVector_d (dim, first, last)
    {}
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

#ifndef VECTOR_WRAPPER
  Self operator*(const RT& c) const
    { return c * (const RVector_d& )(*this) ; }
  Self operator*(const Quotient<RT>& q) const
    {
      return (q.numerator() * (const RVector_d& )(*this)) /
              q.denominator();
    }
  Self operator/(const Quotient<RT>& q) const
    {
      return (q.denominator() * (const RVector_d& )(*this)) /
              q.numerator();
    }
#endif // VECTOR_WRAPPER

  Self operator/(const RT& c) const
    { return (const RVector_d& )(*this) / c; }
  CGAL::Direction_d<R> direction() const
    { return RVector_d::direction(); }
  Self transform(const CGAL::Aff_transformation_d<R>& t) const
    { return RVector_d::transform(t); }

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

#ifdef VECTOR_WRAPPER
template <class T, class R>
_Vector_d_rft_wrapper<R>
multiply(const Quotient<T>& q,
              const CGAL::Vector_d<R>& w,
              const Quotient_tag&)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return _Vector_d_rft_wrapper<R>(
                 CGAL::Vector_d<R>((q.numerator() * (const RVector_d& )(w))
                                  / q.denominator()));
}

template < class R >
_Vector_d_rft_wrapper<R>
multiply(const CGAL::Vector_d<R>& v,
              const CGAL::Vector_d<R>& w,
              const No_number_tag&)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return _Vector_d_rft_wrapper<R>((const RVector_d& )(v)
                                     * (const RVector_d& )(w));
}

template < class T, class R >
_Vector_d_rft_wrapper<R>
multiply(const T& n,
              const CGAL::Vector_d<R>& w,
              const Number_tag&)
{
  typedef typename  R::Vector_d_base  RVector_d;
  typedef typename  R::RT             RT;
  return _Vector_d_rft_wrapper<R>(
                 CGAL::Vector_d<R>(RT(n) * (const RVector_d& )(w)));
}

template <class T, class R>
_Vector_d_rft_wrapper<R>
operator*(const T& t, const CGAL::Vector_d<R>& w)
{
  return multiply(t, w, number_type_tag(t));
}
#endif // VECTOR_WRAPPER

#ifndef NO_OSTREAM_INSERT_VECTOR_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Vector_d<R>& v)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return os << (const RVector_d& )v;
}
#endif // NO_OSTREAM_INSERT_VECTOR_D

#ifndef NO_ISTREAM_EXTRACT_VECTOR_D
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Vector_d<R>& p)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return is >> (RVector_d& )p;
}
#endif // NO_ISTREAM_EXTRACT_VECTOR_D


template<class R>
inline
CGAL::Vector_d<R>
cross_product(const CGAL::Vector_d<R>& v,const CGAL::Vector_d<R>& w)
{
  typedef typename  R::Vector_d_base  RVector_d;
  return cross_product((const RVector_d& )v,(const RVector_d& )w);
}

CGAL_END_NAMESPACE


#endif // CGAL_VECTOR_D_H
