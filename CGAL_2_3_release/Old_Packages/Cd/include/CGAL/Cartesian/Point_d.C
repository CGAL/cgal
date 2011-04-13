// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Point_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_D_C
#define CGAL_CARTESIAN_POINT_D_C

#include <CGAL/Origin.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian/d_utils.h>

#ifndef CGAL_CTAG
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
PointCd<R CGAL_CTAG>::PointCd(int dim)
{
  PTR = new _d_tuple<FT>(dim);
}

template < class R >
PointCd<R CGAL_CTAG>::
PointCd(const PointCd<R CGAL_CTAG> &p)
  : Handle(p)
{}

template < class R >
inline
PointCd<R CGAL_CTAG>::
PointCd(const typename PointCd<R CGAL_CTAG>::Vector_d &v)
  : Handle((const Handle&)v)
{
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointCd<R CGAL_CTAG>::
PointCd(int dim, const Origin &)
{
  CGAL_kernel_precondition (dim>=0);
  PTR = new _d_tuple<FT>(dim);
  std::fill(begin(), end(), FT(0));
}

template < class R >
inline
PointCd<R CGAL_CTAG>::~PointCd()
{}

template < class R >
inline
PointCd<R CGAL_CTAG> &
PointCd<R CGAL_CTAG>::operator=(const PointCd<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
inline
bool
PointCd<R CGAL_CTAG>::operator==(const PointCd<R CGAL_CTAG>& p) const
{
  if (dimension() != p.dimension()) return false;
  if (ptr() == p.ptr()) return true; // identical
  return std::equal(begin(),end(),p.begin());
}

template < class R >
inline
bool 
PointCd<R CGAL_CTAG>::operator!=(const PointCd<R CGAL_CTAG>& p) const
{
  return !(*this == p);
}

template < class R >
inline
bool
PointCd<R CGAL_CTAG>::operator==(const Origin&) const
{
  const_iterator non_zero = find_if(begin(),end(),
                                    bind1st(not_equal_to<FT>(),FT(0)));
  return non_zero == end();
}

template < class R >
inline
bool 
PointCd<R CGAL_CTAG>::operator!=(const Origin&) const
{
  return !(*this == p);
}

template < class R >
inline
long PointCd<R CGAL_CTAG>::id() const
{
  return (long)PTR;
}

template < class R >
CGAL_KERNEL_INLINE
typename PointCd<R CGAL_CTAG>::FT  
PointCd<R CGAL_CTAG>::cartesian(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<dimension()) );
    return *(begin()+i);
}

template < class R >
inline
typename PointCd<R CGAL_CTAG>::FT
PointCd<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition ( (i>=0) && (i<=dimension()) );
  return (i==dimension()) ? FT(1) : cartesian(i);
}

template < class R >
inline
typename PointCd<R CGAL_CTAG>::FT  
PointCd<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
PointCd<R CGAL_CTAG>
PointCd<R CGAL_CTAG>::
transform(const PointCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return t.transform(*this);
}

#ifndef NO_OSTREAM_INSERT_POINTCD
template < class R >
std::ostream& 
operator<<(std::ostream& os, const PointCd<R CGAL_CTAG> &p)
{
  typedef typename R::FT FT;
  print_d<FT> prt(&os);
  if (os.iword(IO::mode)==IO::PRETTY) os << "PointCd(";
  prt(p.dimension());
  if (os.iword(IO::mode)==IO::PRETTY) { os << ", ("; prt.reset(); }
  std::for_each(p.begin(),p.end(),prt);
  if (os.iword(IO::mode)==IO::PRETTY) os << "))";
  return os;
}
#endif // NO_OSTREAM_INSERT_POINTCD

#ifndef NO_ISTREAM_EXTRACT_POINTCD
template < class R >
std::istream& 
operator>>(std::istream& is, PointCd<R CGAL_CTAG> &p)
{
    int d=0, i;
    typedef typename PointCd<R CGAL_CTAG>::FT FT;
    FT *e=0;
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
    if (is)
        p = PointCd<R CGAL_CTAG>(d, e, e+d);
    delete[] e;
    return is;
}

#endif // NO_ISTREAM_EXTRACT_POINTCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_POINT_D_C
