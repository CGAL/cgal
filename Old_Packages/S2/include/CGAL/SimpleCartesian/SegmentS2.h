// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/SegmentS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_SEGMENTS2_H
#define CGAL_SEGMENTS2_H

#include <CGAL/SimpleCartesian/LineS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class SegmentS2
{
public:
                  SegmentS2();
                  SegmentS2(const PointS2<FT>& sp,
                            const PointS2<FT>& ep);

  bool            is_horizontal() const;
  bool            is_vertical() const;
  bool            has_on(const PointS2<FT>& p) const;
  bool            collinear_has_on(const PointS2<FT>& p) const;

  bool            operator==(const SegmentS2<FT>& s) const;
  bool            operator!=(const SegmentS2<FT>& s) const;
  int             id() const;

  const PointS2<FT>&    start() const;
  const PointS2<FT>&    end() const;

  const PointS2<FT>&    source() const;
  const PointS2<FT>&    target() const;

  PointS2<FT>     min() const;
  PointS2<FT>     max() const;
  PointS2<FT>     vertex(int i) const;
  PointS2<FT>     point(int i) const;
  PointS2<FT>     operator[](int i) const;

  FT              squared_length() const;

  DirectionS2<FT> direction() const;
  LineS2<FT>      supporting_line() const;
  SegmentS2<FT>   opposite() const;
  SegmentS2<FT>   transform(const Aff_transformationS2<FT>& t) const;

  bool            is_degenerate() const;
  Bbox_2          bbox() const;

// private:
  PointS2<FT>     s;
  PointS2<FT>     t;
};

template < class FT >
CGAL_KERNEL_CTOR_INLINE
SegmentS2<FT>::SegmentS2() {}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
SegmentS2<FT>::SegmentS2(const PointS2<FT>& sp,
                         const PointS2<FT>& ep)
 : s(sp), t(ep) {}


template < class FT >
inline
bool  
SegmentS2<FT>::operator==(const SegmentS2<FT>& s) const
{ return ( (source() == s.source()) && (target() == s.target()) ); }

template < class FT >
inline
bool  
SegmentS2<FT>::operator!=(const SegmentS2<FT>& s) const
{ return !(*this == s); }

template < class FT >
inline
const PointS2<FT>&  
SegmentS2<FT>::start() const
{ return s; }

template < class FT >
inline
const PointS2<FT>&  
SegmentS2<FT>::end() const
{ return t; }


template < class FT >
inline
const PointS2<FT>&  
SegmentS2<FT>::source() const
{ return s; }

template < class FT >
inline
const PointS2<FT>&  
SegmentS2<FT>::target() const
{ return t; }


template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT>  
SegmentS2<FT>::min() const
{
  return (lexicographically_xy_smaller(source(),target())) ? source()
                                                           : target();
}

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT>  
SegmentS2<FT>::max() const
{
  return (lexicographically_xy_smaller(source(),target())) ? target()
                                                           : source();
}

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT> 
SegmentS2<FT>::vertex(int i) const
{ return (i%2 ==0) ? source() : target(); }

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT> 
SegmentS2<FT>::point(int i) const
{ return (i%2 ==0) ? source() : target(); }

template < class FT >
inline
PointS2<FT> 
SegmentS2<FT>::operator[](int i) const
{ return vertex(i); }

template < class FT >
CGAL_KERNEL_INLINE
FT 
SegmentS2<FT>::squared_length() const
{ return squared_distance(source(), target()); }

template < class FT >
CGAL_KERNEL_INLINE
DirectionS2<FT> 
SegmentS2<FT>::direction() const
{ return DirectionS2<FT>( target() - source() ); }

template < class FT >
inline
LineS2<FT> 
SegmentS2<FT>::supporting_line() const
{ return LineS2<FT>(*this); }

template < class FT >
inline
SegmentS2<FT> 
SegmentS2<FT>::opposite() const
{ return SegmentS2<FT>(target(), source()); }

template < class FT >
inline
SegmentS2<FT> 
SegmentS2<FT>::transform( const Aff_transformationS2<FT>& t) const
{ return SegmentS2<FT>(t.transform(source()), t.transform(target())); }

template < class FT >
CGAL_KERNEL_INLINE
Bbox_2 
SegmentS2<FT>::bbox() const
{ return source().bbox() + target().bbox(); }

template < class FT >
inline
bool  
SegmentS2<FT>::is_degenerate() const
{ return (source() == target()); }


#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTS2
template < class FT >
std::ostream& operator<<(std::ostream &os, const SegmentS2<FT> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentS2(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTS2
template < class FT >
std::istream& operator>>(std::istream &is, SegmentS2<FT> &s)
{
    PointS2<FT> p, q;

    is >> p >> q;

    s = SegmentS2<FT>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTS2

template < class FT >
CGAL_KERNEL_INLINE
bool 
SegmentS2<FT>::is_horizontal() const
{ return source().y() == target().y(); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
SegmentS2<FT>::is_vertical() const
{ return source().x() == target().x(); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
SegmentS2<FT>::has_on(const PointS2<FT>& p) const
{
  return(( p == source() )
         || ( p == target() )
         || ( collinear(source(), p, target())
             && ( DirectionS2<FT>(p - source())
                  !=
                  DirectionS2<FT>(p - target()))
             )
         );
}


template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool 
SegmentS2<FT>::collinear_has_on(const PointS2<FT>& p) const
{
    CGAL_kernel_exactness_precondition( collinear(source(), p, target()) );
    if (CGAL_NTS abs(target().x()-source().x())
      > CGAL_NTS abs(target().y()-source().y())) {
        if (p.x() < source().x())
            return (p.x() >= target().x());
        if (p.x() <= target().x())
            return true;
        return (p.x() == source().x());
    } else {
        if (p.y() < source().y())
            return (p.y() >= target().y());
        if (p.y() <= target().y())
            return true;
        return (p.y() == source().y());
    }
}


CGAL_END_NAMESPACE

#endif
