// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : PointSd.h
// package       : _d
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//                 Bernd Gaertner
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINTSD_H
#define CGAL_POINTSD_H

#ifndef D_TUPLE_H
#include <CGAL/d_tuple.h>
#endif // D_TUPLE_H

CGAL_BEGIN_NAMESPACE

template < class FT >
class DASd;

template < class FT >
class PointSd : public Handle
{
    friend class DASd<FT>;

    public:
    PointSd ();
    PointSd (int dim, const Origin&);
    PointSd (const PointSd<FT> &p);

    template <class InputIterator>
    PointSd (int dim, InputIterator first, InputIterator last)
    {
      CGAL_kernel_precondition( dim >= 0);
      PTR = new _d_tuple<FT>(dim);
      FT* o;
      FT* e = ptr()->e;
      InputIterator i;
      for ( i=first, o=e; ( i < last)&&( o < e+dim ); *(o++) = *(i++) ) {};
      CGAL_kernel_precondition( o == e+dim );
      if ( i < last )
      {
        FT h = *(i++);
        CGAL_kernel_precondition( !CGAL_NTS is_zero (h) );
        CGAL_kernel_precondition( i == last );
        for ( o=e; o < e+dim; *(o++) /= h ) ;
        // if ( h != FT(1) ) { for ( o=e; o < e+dim; *(o++) /= h ) {}; }
      }
    }


    ~PointSd();

    PointSd<FT> &operator=(const PointSd<FT> &p);

      bool operator==(const PointSd<FT> &p) const;
      bool operator!=(const PointSd<FT> &p) const;

    int id() const;

    FT homogeneous (int i) const;
    FT cartesian (int i) const;
    FT operator[] (int i) const;
    const FT* begin() const;
    const FT* end() const;

    int dimension () const;

    private:
     const _d_tuple<FT>* ptr() const;
};

template < class FT >
CGAL_KERNEL_INLINE
const _d_tuple<FT>* PointSd<FT>::ptr() const
{
  return (_d_tuple<FT>*)PTR;
}
CGAL_END_NAMESPACE


#ifndef CGAL_ORIGIN_H
#include <CGAL/Origin.h>
#endif // CGAL_ORIGIN_H
#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif // CGAL_NUMBER_UTILS_H

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointSd<FT>::PointSd()
{
  PTR = new _d_tuple<FT>();
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointSd<FT>::PointSd(int dim, const Origin &)
{
  CGAL_kernel_precondition (dim>=0);
  PTR = new _d_tuple<FT>(dim);
  FT* e = ptr()->e;
  FT* i;
  for (i=e; i<e+dim; *(i++)=FT(0));  // XXX
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointSd<FT>::PointSd(const PointSd<FT> &p)
  : Handle((Handle&)p)
{}


template < class FT >
inline
PointSd<FT>::~PointSd()
{}

template < class FT >
inline
PointSd<FT> &PointSd<FT>::operator=(const PointSd<FT> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class FT >
inline
bool PointSd<FT>::operator==(const PointSd<FT>& p) const
{
  int d = dimension();
  if (d != p.dimension()) return false;
  FT* e = ptr()->e;
  FT* ep = p.ptr()->e;
  FT* i = e;
  FT* j = ep;
  for ( ; i < e+d; ++i, ++j)
      if (*i != *j) return false;
  return true;
}

template < class FT >
inline
bool
PointSd<FT>::operator!=(const PointSd<FT>& p) const
{ return !(*this == p); }

template < class FT >
inline
int PointSd<FT>::id() const
{ return (int)PTR; }

template < class FT >
inline
FT
PointSd<FT>::homogeneous(int i) const
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
PointSd<FT>::cartesian(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<dimension()) );
    return ptr()->e[i];
}

template < class FT >
inline
FT
PointSd<FT>::operator[](int i) const
{ return cartesian(i); }


template < class FT >
inline
int
PointSd<FT>::dimension() const
{ return ptr()->d; }

template < class FT >
inline
const FT*
PointSd<FT>::begin() const
{ return ptr()->e; }

template < class FT >
inline
const FT* PointSd<FT>::end() const
{ return ptr()->e + dimension(); }

#ifndef NO_OSTREAM_INSERT_POINTSD
template < class FT >
std::ostream&
operator<<(std::ostream& os, const PointSd<FT> &p)
{
    int d = p.dimension();
    int i;
    switch(os.iword(IO::mode))
    {
      case IO::ASCII :
        os << d << ' ';
        for (i=0; i<d-1; ++i)
        { os << p.cartesian(i) << ' '; }
        os << p.cartesian(d-1);
        return os;
      case IO::BINARY :
        write(os, d);
        for (i=0; i<d; ++i)
           write(os, p.cartesian(i));
        return os;
      default:
        os << "PointSd(";
        os << d << ", (";
        for (i=0; i<d-1; ++i)
            os << p.cartesian(i) << ", ";
        return os << p.cartesian(d-1) << "))";
    }
}
#endif // NO_OSTREAM_INSERT_POINTSD

#ifndef NO_ISTREAM_EXTRACT_POINTSD
template < class FT >
std::istream&
operator>>(std::istream& is, PointSd<FT> &p)
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
    p = PointSd<FT>( d, e, e+d);
    delete[] e;
    return is;
}

#endif // NO_ISTREAM_EXTRACT_POINTSD
CGAL_END_NAMESPACE


#endif // CGAL_POINTSD_H
