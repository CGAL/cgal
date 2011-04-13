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
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/IO/Pm_Postscript_file_stream.h
// package       : pm (5.45)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_IO_PM_POSTSCRIPT_FILE_STREAM_H
#define CGAL_IO_PM_POSTSCRIPT_FILE_STREAM_H

#ifndef CGAL_LEDA_WINDOW_H
#include <CGAL/IO/leda_window.h>
#endif

#ifndef CGAL_POSTSCRIPT_FILE_STREAM_H
#include <CGAL/IO/Postscript_file_stream.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif
  
/*  #ifndef CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H
  #include <CGAL/IO/Pm_bounding_box_base_Window_stream.h>
  #endif
*/

#ifndef CGAL_IO_PM_DRAWER_H
#include <CGAL/IO/Pm_drawer.h>
#endif

#ifndef CGAL_IO_DRAW_PM_H
#include  <CGAL/IO/draw_pm.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Dcel,class Traits>
Postscript_file_stream& operator << (Postscript_file_stream & ps, Planar_map_2<Dcel,Traits> &pm)
{
  Pm_drawer< Planar_map_2<Dcel,Traits> , Postscript_file_stream>  drawer(ps);
  
  draw_pm(pm, drawer, ps);

  return ps;
}  


/*template <class Dcel,class Traits>
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
        } */ 


CGAL_END_NAMESPACE

#endif







