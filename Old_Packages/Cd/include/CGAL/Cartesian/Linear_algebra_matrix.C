// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_C
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_C

#include <CGAL/Cartesian/d_utils.h>
#include <CGAL/Cartesian/Linear_algebra_matrix.h>
#include <iterator>
#include <algorithm>
#include <numeric>

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class LA >
inline
void
LA_matrixCd<LA>::
new_rep(int n, int m)
{
  _rdim = n;
  _cdim = m;
  PTR = new _d_tuple<Vector *>(n*m);
};

template < class LA >
inline
LA_matrixCd<LA>::
LA_matrixCd(int n)
{
  new_rep(n,n);
}

template < class LA >
inline
LA_matrixCd<LA>::
LA_matrixCd(const std::pair<int,int> &d)
{
  new_rep(d.first,d.second);
}

template < class LA >
inline
LA_matrixCd<LA>::
LA_matrixCd(int n, int m)
{
  new_rep(n,m);
}

template < class LA >
inline
LA_matrixCd<LA>::
LA_matrixCd(const LA_matrixCd<LA> &v)
  : Handle(v)
{}

template < class LA >
inline
LA_matrixCd<LA>::~LA_matrixCd()
{
  iterator i;
  for (i=begin(); i!= end(); ++i)
    delete *i;
}

template < class LA >
inline
LA_matrixCd<LA> &
LA_matrixCd<LA>::operator=(const LA_matrixCd<LA> &v)
{
  Handle::operator=((const Handle &)v);
  return *this;
}

template < class LA >
inline
bool
LA_matrixCd<LA>::operator==(const LA_matrixCd<LA> &M) const
{
  if (dimension() != M.dimension()) return false;
  if (ptr() == M.ptr()) return true; // identical
  return std::equal(begin(),end(),M.begin());
}

template < class LA >
inline
bool
LA_matrixCd<LA>::operator!=(const LA_matrixCd<LA> &v) const
{
  return !(*this == v);
}

template < class LA >
inline
const typename LA_matrixCd<LA>::Vector &
LA_matrixCd<LA>::row(int i) const
{
  return *(*(begin()+i));
}

template < class LA >
inline
std::pair<int,int>
LA_matrixCd<LA>::operator[](int i) const
{
  return std::make_pair(row_dimension(),column_dimension());
}

template < class LA >
inline
const typename LA_matrixCd<LA>::Vector &
LA_matrixCd<LA>::operator[](int i) const
{
  return row(i);
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::operator+(const LA_matrixCd<LA> &M) const
{
  CGAL_kernel_precondition( dimension() == M.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::plus<FT>());
  return w;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::operator-(const LA_matrixCd<LA> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::minus<FT>());
  return w;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::operator-() const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::negate<FT>());
  return v;
}

template < class LA >
inline
typename LA_matrixCd<LA>::FT
LA_matrixCd<LA>::operator*(const LA_matrixCd<LA> &v) const
{
  CGAL_kernel_precondition( column_dimension() == v.row_dimension() );
  return ;
}

typename LA_matrixCd<LA>::FT
LA_matrixCd<LA>::operator*(const Vector &v) const
{
  CGAL_kernel_precondition( row_dimension() == v.dimension() );
  return ;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::
operator*(const typename LA_matrixCd<LA>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(multiplies<FT>(),c));
  return v;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::
operator/(const typename LA_matrixCd<LA>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::divides<FT>(),c));
  return v;
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LA_VECTORCD
template < class LA >
std::ostream &
operator<<(std::ostream &os, const LA_matrixCd<LA> &v)
{
  typedef typename LA_matrixCd<LA>::FT FT;
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
template < class LA >
std::istream &
operator>>(std::istream &is, LA_matrixCd<LA> &v)
{
  typedef typename LA_matrixCd<LA>::FT FT;
  int cdim, rdim, dim;
  FT* w;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      is >> rdim >> cdim;
      dim = cdrim*rdim;
      w = new FT[dim];
      std::copy_n(std::istream_iterator<FT>(is),dim, w);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  v = LA_matrixCd<LA>(rdim,cdim,w,w+dim);
  return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LA_VECTORCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_VECTOR_D_C
