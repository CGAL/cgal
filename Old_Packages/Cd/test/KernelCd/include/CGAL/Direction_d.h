#ifndef CGAL_DIRECTION_D_H
#define CGAL_DIRECTION_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Vector_d.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Direction_d : public _R::Direction_d_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef CGAL::Direction_3<R>              Self;
  typedef typename R::Direction_d_base      RDirection_d;
  typedef typename R::Vector_d              Vector_d;
  typedef const FT*                         const_iterator;

  Direction_d()
    {}
  Direction_d(const Self& v) : RDirection_d( (const RDirection_d& )v )
    {}
  Direction_d(const RDirection_d&  v) : RDirection_d(v)
    {}
  template <class InputIterator>
  Direction_d (int dim, InputIterator first, InputIterator last)
      : RDirection_d (dim, first, last)
    {}
  Self& operator=(const Self& v)
    { RDirection_d::operator=(v); return *this; }
  bool operator==(const Self& v) const
    { return RDirection_d::operator==(v); }
  bool operator!=(const Self& v) const
    { return !(*this == v); }
  bool operator==(const Null_vector& v) const
    { return RDirection_d::operator==(v); }
  bool operator!=(const Null_vector& v) const
    { return !(*this == v); }
  FT delta(int i) const
    { return RDirection_d::delta(i); }
  FT operator[](int i) const
    { return delta(i); }
  int dimension() const
    { return RDirection_d::dimension(); }
  Vector_d to_vector() const
    { return RDirection_d::to_vector(); }
  Self transform(const CGAL::Aff_transformation_d<R>& t) const
    { return RDirection_d::transform(t); }
};

#ifndef NO_OSTREAM_INSERT_DIRECTION_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Direction_d<R>& v)
{
  typedef typename  R::Direction_d_base  RDirection_d;
  return os << (const RDirection_d& )v;
}
#endif // NO_OSTREAM_INSERT_DIRECTION_D

#ifndef NO_ISTREAM_EXTRACT_DIRECTION_D
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Direction_d<R>& p)
{
  typedef typename  R::Direction_d_base  RDirection_d;
  return is >> (RDirection_d& )p;
}
#endif // NO_ISTREAM_EXTRACT_DIRECTION_D


template<class R>
inline
CGAL::Direction_d<R>
cross_product(const CGAL::Direction_d<R>& v,const CGAL::Direction_d<R>& w)
{
  typedef typename  R::Direction_d_base  RDirection_d;
  return cross_product((const RDirection_d& )v,(const RDirection_d& )w);
}

CGAL_END_NAMESPACE


#endif // CGAL_DIRECTION_D_H
