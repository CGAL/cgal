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
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_SEGMENTH3_H
#define CGAL_SEGMENTH3_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Segment_repH3 : public Ref_counted
{
public:
   typedef typename R::Kernel_base::Point_3  Point_3;

   Segment_repH3() {}
   Segment_repH3(const Point_3& sp, const Point_3& ep)
     : startpoint(sp), endpoint(ep) {}

   const Point_3 & start() const { return startpoint; }
   const Point_3 & end()   const { return endpoint; }

private:
   Point_3  startpoint;
   Point_3  endpoint;
};

template < class R >
class Simple_Segment_repH3
{
public:
   typedef typename R::Kernel_base::Point_3  Point_3;

   Simple_Segment_repH3() {}
   Simple_Segment_repH3(const Point_3& sp, const Point_3& ep)
     : startpoint(sp), endpoint(ep) {}

   const Point_3 & start() const { return startpoint; }
   const Point_3 & end()   const { return endpoint; }

private:
   Point_3  startpoint;
   Point_3  endpoint;
};

template < class R_ >
class SegmentH3
  : public R_::Segment_handle_3
{
  typedef typename R_::RT   RT;
  typedef typename R_::FT   FT;
  typedef typename R_::Kernel_base::Point_3              Point_3;
  typedef typename R_::Kernel_base::Direction_3          Direction_3;
  typedef typename R_::Kernel_base::Line_3               Line_3;
  typedef typename R_::Kernel_base::Aff_transformation_3 Aff_transformation_3;

  typedef typename R_::Segment_handle_3            Segment_handle_3_;
  typedef typename Segment_handle_3_::element_type Segment_ref_3;

public:
  typedef R_               R;
 
  SegmentH3()
    : Segment_handle_3_(Segment_ref_3()) {}

  SegmentH3( const Point_3& sp, const Point_3& ep)
    : Segment_handle_3_(Segment_ref_3(sp, ep)) {}

  const Point_3 &  source() const;
  const Point_3 &  target() const;

  const Point_3 &  start() const;
  const Point_3 &  end() const;
  const Point_3 &  min() const;
  const Point_3 &  max() const;
  const Point_3 &  vertex(int i) const;
  const Point_3 &  point(int i) const;
  const Point_3 &  operator[](int i) const;

  FT                squared_length() const;
  Direction_3   direction() const;
  Line_3        supporting_line() const;
  SegmentH3<R>  opposite() const;
  SegmentH3<R>  transform( const Aff_transformation_3& t) const;
  Bbox_3            bbox() const;
  bool              has_on(const Point_3& p) const;
  bool              collinear_has_on(const Point_3& p) const;
  bool              is_degenerate() const;

  bool              operator==(const SegmentH3<R>& s) const;
  bool              operator!=(const SegmentH3<R>& s) const;
};


template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::source() const
{ return Ptr()->start(); }

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::target() const
{ return Ptr()->end(); }

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::start() const
{ return Ptr()->start(); }

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::end() const
{ return Ptr()->end(); }

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::min() const
{
  return
  lexicographically_xyz_smaller(target(),source()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::max() const
{
  return lexicographically_xyz_smaller_or_equal(source(),target()) ?
                                                         target() : source();
}

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::vertex(int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::point(int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }

template < class R >
inline
const typename SegmentH3<R>::Point_3 &
SegmentH3<R>::operator[](int i) const
{ return ( i%2 == 0 ) ? start() : end() ; }


template < class R >
CGAL_KERNEL_INLINE
typename SegmentH3<R>::FT
SegmentH3<R>::squared_length() const
{
  return  (Ptr()->end() - Ptr()->start()) *
          (Ptr()->end() - Ptr()->start())   ;
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentH3<R>::Direction_3
SegmentH3<R>::direction() const
{ return Direction_3( Ptr()->end() - Ptr()->start() ); }

template < class R >
CGAL_KERNEL_INLINE
typename SegmentH3<R>::Line_3
SegmentH3<R>::supporting_line() const
{ return Line_3(Ptr()->start(), Ptr()->end()); }

template < class R >
CGAL_KERNEL_INLINE
SegmentH3<R>
SegmentH3<R>::opposite() const
{ return SegmentH3<R>(Ptr()->end(), Ptr()->start()); }

template < class R >
CGAL_KERNEL_INLINE
SegmentH3<R>
SegmentH3<R>::
transform( const typename SegmentH3<R>::Aff_transformation_3& t) const
{
  return SegmentH3<R>(t.transform(Ptr()->start()),
                      t.transform(Ptr()->end())   );
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
  typename SegmentH3<R>::Point_3 p, q;
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
SegmentH3<R>::has_on(const typename SegmentH3<R>::Point_3 &p) const
{
  return( ( p == start() )
       || ( p == end() )
       || (  ( collinear(p,source(),target() )
           &&( Direction_3( p - Ptr()->start())
               !=
               Direction_3( p - Ptr()->end()))
             )
          )
       );
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentH3<R>::collinear_has_on(const typename SegmentH3<R>::Point_3 &p) const
{
  return( ( p == start() )
       || ( p == end() )
       || ( Direction_3( p - Ptr()->start())
            !=
            Direction_3( p - Ptr()->end())
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
