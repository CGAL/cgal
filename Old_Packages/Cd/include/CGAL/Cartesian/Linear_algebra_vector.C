// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_C
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_C

#include <CGAL/Cartesian/d_utils.h>
#include <iterator>
#include <algorithm>
#include <numeric>

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
LA_vectorCd<R CGAL_CTAG>::
LA_vectorCd(int d)
{
  PTR = new _d_tuple<FT>(d);
}

template < class R >
LA_vectorCd<R CGAL_CTAG>::
LA_vectorCd(const LA_vectorCd<R CGAL_CTAG> &v)
  : Handle(v)
{}

template < class R >
LA_vectorCd<R CGAL_CTAG>::~LA_vectorCd()
{}

template < class R >
LA_vectorCd<R CGAL_CTAG> &
LA_vectorCd<R CGAL_CTAG>::operator=(const LA_vectorCd<R CGAL_CTAG> &v)
{
  Handle::operator=((const Handle &)v);
  return *this;
}

template < class R >
bool
LA_vectorCd<R CGAL_CTAG>::operator==(const LA_vectorCd<R CGAL_CTAG> &v) const
{
  if (dimension() != v.dimension()) return false;
  if (ptr() == v.ptr()) return true; // identical
  return std::equal(begin(),end(),v.begin());
}

template < class R >
inline
bool
LA_vectorCd<R CGAL_CTAG>::operator!=(const LA_vectorCd<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}

template < class R >
inline
typename LA_vectorCd<R CGAL_CTAG>::FT
LA_vectorCd<R CGAL_CTAG>::operator[](int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<dimension()) );
  return *(begin()+i);
}

template < class R >
inline
LA_vectorCd<R CGAL_CTAG>
LA_vectorCd<R CGAL_CTAG>::operator+(const LA_vectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::plus<FT>());
  return w;
}

template < class R >
inline
LA_vectorCd<R CGAL_CTAG>
LA_vectorCd<R CGAL_CTAG>::operator-(const LA_vectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::minus<FT>());
  return w;
}

template < class R >
inline
LA_vectorCd<R CGAL_CTAG>
LA_vectorCd<R CGAL_CTAG>::operator-() const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::negate<FT>());
  return v;
}

template < class R >
inline
typename LA_vectorCd<R CGAL_CTAG>::FT
LA_vectorCd<R CGAL_CTAG>::operator*(const LA_vectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  return std::inner_product(begin(),end(),v.begin(),FT(0));
}

template < class R >
inline
LA_vectorCd<R CGAL_CTAG>
LA_vectorCd<R CGAL_CTAG>::
operator*(const typename LA_vectorCd<R CGAL_CTAG>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::multiplies<FT>(),c));
  return v;
}

template < class R >
inline
LA_vectorCd<R CGAL_CTAG>
LA_vectorCd<R CGAL_CTAG>::
operator/(const typename LA_vectorCd<R CGAL_CTAG>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::divides<FT>(),c));
  return v;
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LA_VECTORCD
template < class R >
std::ostream &
operator<<(std::ostream &os, const LA_vectorCd<R CGAL_CTAG> &v)
{
  typedef typename LA_vectorCd<R CGAL_CTAG>::FT FT;
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      os << v.dimension();
      std::for_each(v.begin(),v.end(),print_d<FT>(os));
      break;
    default:
      os << v.dimension() << ", ";
      std::for_each(v.begin(),v.end(),print_d<FT>(os));
  }
  return os;
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_LA_VECTORCD

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LA_VECTORCD
template < class R >
std::istream &
operator>>(std::istream &is, LA_vectorCd<R CGAL_CTAG> &v)
{
  typedef typename LA_vectorCd<R CGAL_CTAG>::FT FT;
  int dim;
  FT* w;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      is >> dim;
      w = new FT[dim];
      std::copy_n(std::istream_iterator<FT>(is),dim, w);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  v = LA_vectorCd<R CGAL_CTAG>(dim,w,w+dim);
  return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LA_VECTORCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_VECTOR_D_C
