#ifndef CGAL_CARTESIAN_POINT_2_C
#define CGAL_CARTESIAN_POINT_2_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#ifndef CGAL_ORIGIN_H
#include <CGAL/Origin.h>
#endif // CGAL_ORIGIN_H
#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif // CGAL_NUMBER_UTILS_H

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_INLINE
const _d_tuple<FT>* PointCd<FT>::ptr() const
{
  return (_d_tuple<FT>*)PTR;
}

template <class R, class InputIterator>
PointCd<R>::
PointCd CGAL_NULL_TMPL_ARGS (int dim, InputIterator first, InputIterator last)
{
  CGAL_kernel_precondition( dim >= 0);
  CGAL_kernel_precondition( (last-first == dim) || (last-first == dim+1) );
  PTR = new _d_tuple<FT>(dim);
  FT* o, e = ptr()->e;
  InputIterator i;
  for (i=first, o=e; (i != last) && (o < e+dim); ++i, ++o ) *o = *i;
  if (i != last) {
    FT h = *(i++);
    CGAL_kernel_precondition( !is_zero (h) );
    CGAL_kernel_precondition( i == last );
    for ( o=e; o < e+dim; o++ ) *o /= h;
  }
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointCd<FT>::PointCd()
{
  PTR = new _d_tuple<FT>();
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointCd<FT>::PointCd(int dim, const Origin &)
{
  CGAL_kernel_precondition (dim>=0);
  PTR = new _d_tuple<FT>(dim);
  FT *e = ptr()->e, *i;
  for (i=e; i<e+dim; *(i++)=FT(0));
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointCd<FT>::PointCd(const PointCd<FT> &p)
  : Handle((Handle&)p)
{}


template < class FT >
inline
PointCd<FT>::~PointCd()
{}

template < class FT >
inline
PointCd<FT> &PointCd<FT>::operator=(const PointCd<FT> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class FT >
inline
bool PointCd<FT>::operator==(const PointCd<FT>& p) const
{
  int d = dimension();
  if (d != p.dimension()) return false;
  FT *e = ptr()->e, *ep = p.ptr()->e, *i, *j;
  for (i=e, j=ep; i<e+d; ++i, ++j)
      if (*i != *j) return false;
  return true;
}

template < class FT >
inline
bool 
PointCd<FT>::operator!=(const PointCd<FT>& p) const
{ return !(*this == p); }

template < class FT >
inline
int PointCd<FT>::id() const
{ return (int)PTR; }

template < class FT >
inline
FT  
PointCd<FT>::homogeneous(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<=dimension()) );
    if (i<dimension())
    return ptr()->e[i];
    else
    return FT(1);
}

template < class FT >
CGAL_KERNEL_INLINE
FT  
PointCd<FT>::cartesian(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<dimension()) );
    return ptr()->e[i];
}

template < class FT >
inline
FT  
PointCd<FT>::operator[](int i) const
{ return cartesian(i); }


template < class FT >
inline
int 
PointCd<FT>::dimension() const
{ return ptr()->d; }

template < class FT >
inline
const FT* 
PointCd<FT>::begin() const
{ return ptr()->e; }

template < class FT >
inline
const FT* PointCd<FT>::end() const
{ return ptr()->e + dimension(); }

#ifndef NO_OSTREAM_INSERT_POINTCD
template < class FT >
std::ostream& 
operator<<(std::ostream& os, const PointCd<FT> &p)
{
    int d = p.dimension(), i;
    switch(os.iword(IO::mode)) 
    {
      case IO::ASCII : 
	os << d << ' ';
        for (i=0; i<d-1; ++i)
          os << p.cartesian(i) << ' ';
        os << p.cartesian(d-1); 
        return os;
      case IO::BINARY :
        write(os, d);
        for (i=0; i<d; ++i)
           write(os, p.cartesian(i));
        return os;
      default:
        os << "PointCd(";
        os << d << ", (";
        for (int i=0; i<d-1; ++i)
            os << p.cartesian(i) << ", ";
        return os << p.cartesian(d-1) << "))";
    }
}
#endif // NO_OSTREAM_INSERT_POINTCD

#ifndef NO_ISTREAM_EXTRACT_POINTCD
template < class FT >
std::istream& 
operator>>(std::istream& is, PointCd<FT> &p)
{
    int d=0, i;
    FT* e=0;
    switch(is.iword(IO::mode)) 
    {
    case IO::ASCII :
      is >> d;
      e = new FT[d];
      for (i=0; i < d; ++i) { is >> e[i]; }
      break;
    case IO::BINARY :
      read(is, d);
      e = new FT[d];
      for (i=0; i<d; ++i) { read(is, e[i]); }
      break;
    default:
      CGAL_kernel_assertion_msg(false,\
                                "Stream must be in ascii or binary mode"); 
      // throw ios_base::failure("Stream must be in ascii or binary mode");
      return is;
    }
    p = PointCd<FT>( d, e, e+d);
    delete[] e;
    return is;
}

#endif // NO_ISTREAM_EXTRACT_POINTCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_POINT_D_C
