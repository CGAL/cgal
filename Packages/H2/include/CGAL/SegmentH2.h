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
// file          : SegmentH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_SEGMENTH2_H
#define CGAL_SEGMENTH2_H

#include <CGAL/PointH2.h>
#include <CGAL/LineH2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Segment_repH2 : public Ref_counted
{
public:
            Segment_repH2();
            Segment_repH2(const PointH2<R>& sp,
                          const PointH2<R>& ep);

    PointH2<R>  start;
    PointH2<R>  end;
};

template < class R_ >
class SegmentH2
  : public R_::Segment_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

            SegmentH2();
            SegmentH2( const PointH2<R>& sp,
                       const PointH2<R>& ep);
            SegmentH2( const RT& sw, const RT& sx, const RT& sy,
                       const RT& ew, const RT& ex, const RT& ey);

    bool    operator==(const SegmentH2<R>& s) const;
    bool    operator!=(const SegmentH2<R>& s) const;

    PointH2<R>  source() const;
    PointH2<R>  target() const;
    PointH2<R>  start() const;
    PointH2<R>  end()   const;
    PointH2<R>  vertex(int i) const;
    PointH2<R>  point(int i) const;
    PointH2<R>  operator[](int i) const;
    PointH2<R>  min()   const;
    PointH2<R>  max()   const;
    PointH2<R>  other_vertex(const PointH2<R>& p) const;

    bool    is_horizontal() const;
    bool    is_vertical() const;
    bool    has_on(const PointH2<R>& p) const;
    bool    collinear_has_on(const PointH2<R>& p) const;
    bool    is_degenerate() const;

    FT      squared_length() const;

    DirectionH2<R> direction() const;
    LineH2<R>      supporting_line() const;
    SegmentH2<R>   opposite() const;
    Bbox_2             bbox() const;

    SegmentH2<R>
            transform( const Aff_transformationH2<R> & t) const;

};



template < class R >
CGAL_KERNEL_CTOR_INLINE
Segment_repH2<R>::Segment_repH2()
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
Segment_repH2<R>::
Segment_repH2(const PointH2<R>& sp, const PointH2<R>& ep)
 : start(sp), end(ep)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentH2<R>::SegmentH2()
 : Handle_for< Segment_repH2<R> >( Segment_repH2<R>() )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentH2<R>::SegmentH2( const PointH2<R>& sp,
                             const PointH2<R>& ep)
 : Handle_for< Segment_repH2<R> >( Segment_repH2<R>(sp,ep) )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentH2<R>::SegmentH2(const RT& sx, const RT& sy, const RT& sw,
                            const RT& ex, const RT& ey, const RT& ew)
 : Handle_for< Segment_repH2<R> >( Segment_repH2<R>(
                                         PointH2<R>(sx,sy,sw),
                                         PointH2<R>(ex,ey,ew) ) )
{}

template < class R >
inline
PointH2<R>
SegmentH2<R>::source() const
{ return ptr->start; }

template < class R >
inline
PointH2<R>
SegmentH2<R>::start() const
{ return ptr->start; }

template < class R >
inline
PointH2<R>
SegmentH2<R>::target() const
{ return ptr->end; }

template < class R >
inline
PointH2<R>
SegmentH2<R>::end() const
{ return ptr->end; }

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
SegmentH2<R>::min() const
{
  return
  lexicographically_xy_smaller_or_equal(start(),end()) ? start() : end();
}

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
SegmentH2<R>::max() const
{
  return
  lexicographically_xy_smaller_or_equal(start(),end()) ? end() : start();
}

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
SegmentH2<R>::other_vertex(const PointH2<R>& p) const
{
  CGAL_kernel_precondition( (p == end()) || (p == start()) );
  return ( p == start() ) ? end() : start() ;
}

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
SegmentH2<R>::vertex(int i) const
{
  switch (i%2)
  {
    case 0:  return ptr->start;
    case 1:  return ptr->end;
  };
  return PointH2<R>(); // otherwise the SGI compiler complains
}

template < class R >
inline
PointH2<R>
SegmentH2<R>::point(int i) const
{ return vertex(i); }

template < class R >
inline
PointH2<R>
SegmentH2<R>::operator[](int i) const
{ return vertex(i); }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
SegmentH2<R>::squared_length() const
{ return  (ptr->end - ptr->start) * (ptr->end - ptr->start); }

template < class R >
CGAL_KERNEL_INLINE
DirectionH2<R>
SegmentH2<R>::direction() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return DirectionH2<R>( ptr->end - ptr->start );
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>
SegmentH2<R>::supporting_line() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return LineH2<R>(ptr->start, ptr->end);
}

template < class R >
CGAL_KERNEL_INLINE
SegmentH2<R>
SegmentH2<R>::opposite() const
{ return SegmentH2<R>(ptr->end, ptr->start); }

template < class R >
CGAL_KERNEL_INLINE
SegmentH2<R>
SegmentH2<R>::
transform(const Aff_transformationH2<R>& t) const
{
  return SegmentH2<R>(t.transform(ptr->start),
                               t.transform(ptr->end)   );
}

template < class R >
inline
Bbox_2
SegmentH2<R>::bbox() const
{ return start().bbox() + end().bbox(); }

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentH2<R> &s)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentH2(" << s.source() <<  ", " << s.target() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTH2
template < class R >
std::istream &
operator>>(std::istream &is, SegmentH2<R> &s)
{
  PointH2<R> p, q;
  is >> p >> q;
  s = SegmentH2<R>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTH2

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::is_horizontal() const
{
  return (    ptr->start.hy() * ptr->end.hw()
           == ptr->end.hy() * ptr->start.hw() );
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::is_vertical() const
{
  return (    ptr->start.hx() * ptr->end.hw()
           == ptr->end.hx() * ptr->start.hw() );
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::is_degenerate() const
{ return (start() == end()); }
template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::has_on(const PointH2<R>& p) const
{
  if ( collinear(ptr->start, p, ptr->end ) )
  {
      return collinear_has_on(p);
  }
  else
  {
      return false;
  }
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::collinear_has_on(const PointH2<R>& p) const
{
  return (  lexicographically_xy_smaller_or_equal(p, max() )
         && lexicographically_xy_smaller_or_equal(min(),p) );
}
template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH2<R>::operator==(const SegmentH2<R>& s) const
{
  return (  (start() == s.start() )
          &&(end()   == s.end()   ) );
}

template < class R >
inline
bool
SegmentH2<R>::operator!=(const SegmentH2<R>& s) const
{ return ( !operator==(s) ); }

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENTH2_H
