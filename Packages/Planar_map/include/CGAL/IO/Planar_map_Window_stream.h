// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>

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


