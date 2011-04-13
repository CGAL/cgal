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
// file          : include/CGAL/IO/Arr_Postscript_file_stream.h
// package       : Arrangement (1.82)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_IO_ARR_POSTSCRIPT_FILE_STREAM_H
#define CGAL_IO_ARR_POSTSCRIPT_FILE_STREAM_H

#ifndef CGAL_POSTSCRIPT_FILE_STREAM_H
#include <CGAL/IO/Postscript_file_stream.h>
#endif

#ifndef CGAL_ARRANGEMENT_2_H
#include <CGAL/Arrangement_2.h>
#endif

/*
#ifndef CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H
#include <CGAL/IO/Pm_bounding_box_base_Window_stream.h>
#endif
*/

#ifndef CGAL_IO_FILE_DRAWER_H
#include <CGAL/IO/Pm_drawer.h>
#endif

#ifndef CGAL_IO_DRAW_PM_H
#include <CGAL/IO/draw_pm.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Dcel,class Traits, class Base_node>
Postscript_file_stream& operator << (Postscript_file_stream & ps,  
                                     const Arrangement_2<Dcel, Traits, 
                                     Base_node>& arr)
{

  Pm_drawer<  Arrangement_2<Dcel,Traits, Base_node> , 
    Postscript_file_stream>  drawer(ps);
  
  draw_pm(arr, drawer, ps);

  return ps;
}  

CGAL_END_NAMESPACE

#endif

