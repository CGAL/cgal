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
// file          : include/CGAL/IO/Pm_bounding_box_base_Window_stream.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H
#define CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

#ifndef CGAL_LEDA_WINDOW_H
#include <CGAL/leda_window.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map>
inline Window_stream& operator<<(Window_stream& os,
                          const Pm_bounding_box_base<Planar_map> &b)
{}  

CGAL_END_NAMESPACE

#endif // CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H


