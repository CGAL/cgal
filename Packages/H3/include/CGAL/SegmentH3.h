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
// file          : SegmentH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_SEGMENTH3_H
#define CGAL_SEGMENTH3_H

#include <CGAL/LineH3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Segment_repH3 : public Ref_counted
{
public:
   Segment_repH3() {}
   Segment_repH3(const PointH3<R>& sp, const PointH3<R>& ep)
     : startpoint(sp), endpoint(ep)
   {}

   PointH3<R>  start() const;
   PointH3<R>  end()   const;

private:
   PointH3<R>  startpoint;
   PointH3<R>  endpoint;
};

template < class R_ >
class SegmentH3
  : public R_::Segment_handle_3
{
public:
  typedef R_               R;
  typedef typename R::RT   RT;
  typedef typename R::FT   FT;

  SegmentH3();
  SegmentH3( const PointH3<R>& sp, const PointH3<R>& ep);

  PointH3<R>    source() const;
  PointH3<R>    target() const;

  PointH3<R>    start() const;
  PointH3<R>    end() const;
  PointH3<R>    min() const;
  PointH3<R>    max() const;
  PointH3<R>    vertex(int i) const;
  PointH3<R>    point(int i) const;
  PointH3<R>    operator[](int i) const;

  FT                squared_length() const;
  DirectionH3<R>
                    direction() const;
  LineH3<R>     supporting_line() const;
  SegmentH3<R>  opposite() const;
  SegmentH3<R>  transform( const Aff_transformationH3<R> & t) const;
  Bbox_3            bbox() const;
  bool              has_on(const PointH3<R> p) const;
  bool              collinear_has_on(const PointH3<R> p) const;
  bool              is_degenerate() const;

  bool              operator==(const SegmentH3<R>& s) const;
  bool              operator!=(const SegmentH3<R>& s) const;

};

template < class R >
inline
PointH3<R>
Segment_repH3<R>::start() const
{ return startpoint; }

template < class R >
inline
PointH3<R>
Segment_repH3<R>::end() const
{ return endpoint; }


template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentH3<R>::SegmentH3()
 : Handle_for< Segment_repH3<R> >( Segment_repH3<R>() )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SegmentH3<R>::SegmentH3( const PointH3<R>& sp,
                             const PointH3<R>& ep)
 : Handle_for< Segment_repH3<R> >( Segment_repH3<R>(sp,ep) )
{}
template < class R >
inline
PointH3<R>
SegmentH3<R>::source() const
{ return ptr->start(); }

template < class R >
inline
PointH3<R>
SegmentH3<R>::target() const
{ return ptr->end(); }

template < class R >
inline
PointH3<R>
SegmentH3<R>::start() const
{ return ptr->start(); }

template < class R >
inline
PointH3<R>
SegmentH3<R>::end() const
{ return ptr->end(); }

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
SegmentH3<R>::min() const
{
  return
  lexicographically_xyz_smaller(target(),source()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
SegmentH3<R>::max() const
{
  return
  lexicographically_xyz_smaller_or_equal(source(),target()) ?
                                                         target() : source();
}

template < class R >
inline
PointH3<R>
SegmentH3<R>::vertex(int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }

template < class R >
inline
PointH3<R>
SegmentH3<R>::point(int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }

template < class R >
inline
PointH3<R>
SegmentH3<R>::operator[](int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }


template < class R >
CGAL_KERNEL_INLINE
typename R::FT
SegmentH3<R>::squared_length() const
{
  return  (ptr->end() - ptr->start()) *
          (ptr->end() - ptr->start())   ;
}

template < class R >
CGAL_KERNEL_INLINE
DirectionH3<R>
SegmentH3<R>::direction() const
{ return DirectionH3<R>( ptr->end() - ptr->start() ); }

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
SegmentH3<R>::supporting_line() const
{ return LineH3<R>(ptr->start(), ptr->end()); }

template < class R >
CGAL_KERNEL_INLINE
SegmentH3<R>
SegmentH3<R>::opposite() const
{ return SegmentH3<R>(ptr->end(), ptr->start()); }

template < class R >
CGAL_KERNEL_INLINE
SegmentH3<R>
SegmentH3<R>::
transform( const Aff_transformationH3<R>& t) const
{
  return SegmentH3<R>(t.transform(ptr->start()),
                               t.transform(ptr->end())   );
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_3
SegmentH3<R>::bbox() const
{ return source().bbox() + target().bbox(); }


#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTH3
template < class R >
std::ostream &operator<<(std::ostream &os, const SegmentH3<R> &s)
{
  switch(os.iword(IO::mode))
  {
      case IO::ASCII :
          return os << s.source() << ' ' << s.target();
      case IO::BINARY :
          return os << s.source() << s.target();
      default:
          return os << "SegmentH3(" << s.source()
	            <<  ", " << s.target() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTH3
template < class R >
std::istream &operator>>(std::istream &is, SegmentH3<R> &s)
{
  PointH3<R> p, q;
  is >> p >> q;
  s = SegmentH3<R>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTH3

template < class R >
inline
bool
SegmentH3<R>::is_degenerate() const
{ return  source()==target(); }

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH3<R>::has_on(const PointH3<R> p) const
{
  return( ( p == start() )
       || ( p == end() )
       || (  ( collinear(p,source(),target() )
           &&( DirectionH3<R>( p - ptr->start())
               !=
               DirectionH3<R>( p - ptr->end()))
             )
          )
       );
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH3<R>::collinear_has_on(const PointH3<R> p) const
{
  return( ( p == start() )
       || ( p == end() )
       || ( DirectionH3<R>( p - ptr->start())
            !=
            DirectionH3<R>( p - ptr->end())
          )
        );
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH3<R>::operator==(const SegmentH3<R>& s) const
{
  return ( (start() == s.start() )
         &&(end()   == s.end()   ) );
}

template < class R >
inline
bool
SegmentH3<R>::operator!=(const SegmentH3<R>& s) const
{ return ( !operator==(s) ); }

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENTH3_H
