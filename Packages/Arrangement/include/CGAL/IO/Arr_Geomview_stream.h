// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-62 $
// release_date  : $CGAL_Date: 2001/05/11 $
//
// file          : include/CGAL/IO/Arr_Geomview_stream.h
// package       : Arrangement (1.82)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_IO_ARR_GEOMVIEW_STREAM_H
#define CGAL_IO_ARR_GEOMVIEW_STREAM_H

#ifndef CGAL_ARRANGEMENT_2_H
#include <CGAL/Arrangement_2.h>
#endif

#ifndef CGAL_GEOMVIEW_STREAM_H
#include <CGAL/IO/Geomview_stream.h>
#endif

#ifndef CGAL_IO_FILE_DRAWER_H
#include <CGAL/IO/Pm_drawer.h>
#endif

#ifndef CGAL_IO_DRAW_PM_H
#include <CGAL/IO/draw_pm.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Dcel,class Traits, class Base_node>
Geomview_stream& operator << (Geomview_stream& os, 
                              const Arrangement_2<Dcel,Traits, Base_node>& arr)
{

  Pm_drawer< Arrangement_2<Dcel,Traits, Base_node> , Geomview_stream>  
                                                              drawer(os);
  
  draw_pm(arr, drawer, os);

  return os;
}

CGAL_END_NAMESPACE

#endif




