
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
// file          : include/CGAL/IO/polygonal_2.h
// package       : window
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_IO_POLYGONAL_2_H
#define CGAL_IO_POLYGONAL_2_H

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/stl_extensions.h>
#include <CGAL/ch_value_type.h>
#include <CGAL/circulator.h>

namespace CGAL {

//---------------------- Polygon --------------------------

template <class Stream, class ForwardIterator, class Traits>
void
send_to_stream_as_polygon(Stream& W,
                          ForwardIterator first, ForwardIterator last,
                          const Traits& );


template <class Stream, class ForwardIterator, class R>
inline
void
_send_to_stream_as_polygon(Stream& W,
                           ForwardIterator first, ForwardIterator last,
                           Point_2<R>* )
{ send_to_stream_as_polygon(W, first, last, convex_hull_traits_2<R>() ); }


template <class Stream, class ForwardIterator>
inline
void
send_to_stream_as_polygon(Stream& W,
                          ForwardIterator first, ForwardIterator last)
{ _send_to_stream_as_polygon(W, first, last, ch_value_type(first)); }


template <class Stream, class ForwardIterator, class Traits>
void
send_to_stream_as_polygon(Stream& W,
                          ForwardIterator first, ForwardIterator last,
                          const Traits& )
{
  typedef  typename Traits::Segment_2   Segment2;
  if (first == last) return;
  ForwardIterator it   = first;
  ForwardIterator fifi = CGAL::successor(first);
  while ( fifi != last )
  {
      W << Segment2(*it,*fifi);
      it = fifi++;
  }
  W << Segment2(*it,*first);
  return;
}


//---------------------- Polyline --------------------------


template <class Stream, class ForwardIterator, class Traits>
void
send_to_stream_as_polyline(Stream& W,
                           ForwardIterator first, ForwardIterator last,
                           const Traits& );


template <class Stream, class ForwardIterator, class R>
inline
void
_send_to_stream_as_polyline(Stream& W,
                            ForwardIterator first, ForwardIterator last,
                            Point_2<R>* )
{ send_to_stream_as_polyline(W, first, last, convex_hull_traits_2<R>() ); }


template <class Stream, class ForwardIterator>
inline
void
send_to_stream_as_polyline(Stream& W,
                           ForwardIterator first, ForwardIterator last)
{ _send_to_stream_as_polyline(W, first, last, ch_value_type(first)); }


template <class Stream, class ForwardIterator, class Traits>
void
send_to_stream_as_polyline(Stream& W,
                           ForwardIterator first, ForwardIterator last,
                           const Traits& )
{
  typedef  typename Traits::Segment_2   Segment2;
  if (first == last) return;
  ForwardIterator it   = first;
  ForwardIterator fifi = CGAL::successor(first);
  while ( fifi != last )
  {
      W << Segment2(*it,*fifi);
      it = fifi++;
  }
  return;
}



template <class Stream, class Circulator, class Traits>
void
send_to_stream_as_polygon(Stream& W,
                          Circulator C,
                          const Forward_circulator_tag&,
                          const Traits& );


template <class Stream, class Circulator, class R>
inline
void
_send_to_stream_as_polygon(Stream& W,
                           Circulator C,
                           const Forward_circulator_tag& ct,
                           Point_2<R>* )
{ send_to_stream_as_polygon(W, C, ct, R() ); }


template <class Stream, class Circulator>
inline
void
send_to_stream_as_polygon(Stream& W,
                          Circulator C,
                          const Forward_circulator_tag& ct)
{ _send_to_stream_as_polygon(W, C, ct, ch_value_type(C)); }


template <class Stream, class Circulator, class Traits>
void
send_to_stream_as_polygon(Stream& W,
                          Circulator C,
                          const Forward_circulator_tag&,
                          const Traits& )
{
  typedef  typename Traits::Segment_2   Segment2;
  if ( is_empty_range(C,C) ) return;
  Circulator start = C;
  do
  {
    W << Segment2(*C, *successor(C));
    ++C;
  } while ( C != start );
}


} // namespace CGAL


#endif // CGAL_IO_POLYGONAL_2_H
