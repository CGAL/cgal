// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/IO/Planar_map_Window_stream.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_IO_PLANAR_MAP_WINDOW_STREAM_H
#define CGAL_IO_PLANAR_MAP_WINDOW_STREAM_H

#ifndef CGAL_LEDA_WINDOW_H
#include <CGAL/leda_window.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif
/*
#ifndef CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H
#include <CGAL/IO/Pm_bounding_box_base_Window_stream.h>
#endif
*/
CGAL_BEGIN_NAMESPACE

template <class Dcel,class Traits>
Window_stream& operator<<(Window_stream& os,
                          Planar_map_2<Dcel,Traits> &m)
{
//  os << *m.get_bounding_box();
  Halfedge_iterator it = m.halfedges_begin(), end = m.halfedges_end();

  while(it != end){
    os << it->curve();
    os << it->target()->point();
    os << it->source()->point();
    ++it;++it;
  }
  return os;
}  


template <class Dcel,class Traits>
Window_stream& write(Window_stream& os, Planar_map_2<Dcel,Traits> &m)
{
//  os << *m.get_bounding_box();
  Halfedge_iterator it = m.halfedges_begin(), end = m.halfedges_end();
  const Traits& traits=m.get_traits();
  while(it != end){
	write(os,it->curve(),traits);
    ++it;++it;
  }
  return os;
}  


CGAL_END_NAMESPACE

#endif


