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
// file          : PointHd.h
// package       : _d
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//                 Bernd Gaertner
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINTHD_H
#define CGAL_POINTHD_H

#include <CGAL/d_tuple.h>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT>
class DAHd;

template < class FT, class RT >
class PointHd : public Handle
{
    friend class DAHd<FT,RT>;

    public:
    PointHd ();
    PointHd (int dim, const Origin&);
    PointHd (const PointHd<FT,RT> &p);

    template <class InputIterator>
    PointHd (int dim, InputIterator first, InputIterator last)
    {
      CGAL_kernel_precondition( dim >= 0 );
      PTR = new _d_tuple<RT>( dim, false);
      RT* o;
      RT* e = ptr()->e;
      InputIterator i;
      for (i = first, o = e; ( i < last )&&( o < e+dim ); *(o++) = *(i++) ) {};
      CGAL_kernel_precondition ( o == e+dim );
      if (i < last)
      {
        RT h = *(i++);
        CGAL_kernel_precondition( !CGAL_NTS is_zero (h) );
        CGAL_kernel_precondition( i == last );
        *o = h;
      }
      else
      { *o = RT(1); }
    }

    ~PointHd();

    PointHd<FT,RT>& operator=(const PointHd<FT,RT> &p);

      bool operator==(const PointHd<FT,RT> &p) const;
      bool operator!=(const PointHd<FT,RT> &p) const;

    int id() const;

    RT homogeneous (int i) const;
    FT cartesian (int i) const;
    FT operator[] (int i) const;
    const RT* begin() const;
    const RT* end() const;

    int dimension () const;

    private:
    const _d_tuple<RT>* ptr() const;
};

template < class FT, class RT >
CGAL_KERNEL_INLINE
const _d_tuple<RT>* PointHd<FT,RT>::ptr() const
{
  return (_d_tuple<RT>*)PTR;
}

CGAL_END_NAMESPACE

#include <CGAL/Origin.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointHd<FT,RT>::PointHd()
{ PTR = new _d_tuple<RT>(); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointHd<FT,RT>::PointHd(int dim, const Origin &)
{
  CGAL_kernel_precondition( dim >= 0);
  PTR = new _d_tuple<RT>(dim,false);
  RT* e = ptr()->e;
  FT* i;
  for (i=e; i<e+dim; *(i++)=RT(0));
  *i = RT(1);
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointHd<FT,RT>::PointHd(const PointHd<FT,RT> &p)
  : Handle((Handle&)p)
{}

template < class FT, class RT >
inline
PointHd<FT,RT>::~PointHd()
{}

template < class FT, class RT >
inline
PointHd<FT,RT>&
PointHd<FT,RT>::operator=(const PointHd<FT,RT> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class FT, class RT >
inline
bool
PointHd<FT,RT>::operator==(const PointHd<FT,RT>& p) const
{
  int d = dimension();
  if (d != p.dimension()) return false;
  RT* e = ptr()->e;
  RT* ep = p.ptr()->e;
  RT h = e[d], hp = ep[d];
  RT* i;
  RT* j;
  for ( i=e, j=ep; i<e+d; ++i, ++j)
      if ((*i)*hp != (*j)*h) return false;
  return true;
}

template < class FT, class RT >
inline
bool
PointHd<FT,RT>::operator!=(const PointHd<FT,RT>& p) const
{ return !(*this == p); }

template < class FT, class RT >
inline
int
PointHd<FT,RT>::id() const
{ return (int)PTR; }

template < class FT, class RT >
inline
RT
PointHd<FT,RT>::homogeneous(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<=dimension()) );
    return ptr()->e[i];
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointHd<FT,RT>::cartesian(int i) const
{
    CGAL_kernel_precondition ( (i>=0) && (i<dimension()) );
    return FT(ptr()->e[i])/FT(ptr()->e[dimension()]);
}

template < class FT, class RT >
inline
FT
PointHd<FT,RT>::operator[](int i) const
{ return cartesian(i); }

template < class FT, class RT >
inline
int
PointHd<FT,RT>::dimension() const
{ return ptr()->d; }

template < class FT, class RT >
inline
const RT*
PointHd<FT,RT>::begin() const
{ return ptr()->e; }

template < class FT, class RT >
inline
const RT*
PointHd<FT,RT>::end() const
{ return ptr()->e + dimension() + 1; }

#ifndef CGAL_NO_OSTREAM_INSERT_POINTHD
template < class FT, class RT >
std::ostream&
operator<<(std::ostream& os, const PointHd<FT,RT>& p)
{
    int d = p.dimension();
    int i;
    switch(os.iword(IO::mode))
    {
      case IO::ASCII :
        os << d << ' ';
        for (i=0; i<d; ++i) { os << p.homogeneous(i) << ' '; }
        os << p.homogeneous(d);
        return os;
      case IO::BINARY :
        write(os, d);
        for (i=0; i<d+1; ++i) { write(os, p.homogeneous(i)); }
        return os;
      default:
        os << "PointHd(";
        os << d << ", (";
        for (i=0; i<d; ++i) { os << p.homogeneous(i) << ", "; }
        os << p.homogeneous(d) << "))";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTHD

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTHD
template < class FT, class RT >
std::istream&
operator>>(std::istream& is, PointHd<FT,RT> &p)
{
    int d=0, i;
    RT* e;
    switch(is.iword(IO::mode))
    {
      case IO::ASCII :
        is >> d;
        e = new RT[d+1];
        for (i=0; i<d+1; ++i) { is >> e[i]; }
        break;
      case IO::BINARY :
        read(is, d);
        e = new RT[d+1];
        for (i=0; i<d+1; ++i) { read(is, e[i]); }
        break;
      default:
        CGAL_kernel_assertion_msg(false,\
                                  "Stream must be in ascii or binary mode");
        // throw ios_base::failure("Stream must be in ascii or binary mode");
        return is;
    }
    p = PointHd<FT,RT>( d, e, e+d+1);
    delete[] e;
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTHD

CGAL_END_NAMESPACE

#endif // CGAL_POINTHD_H
