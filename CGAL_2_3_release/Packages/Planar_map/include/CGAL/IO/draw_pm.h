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
// file          : include/CGAL/IO/draw_pm.h
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
#ifndef CGAL_IO_DRAW_PM_H
#define CGAL_IO_DRAW_PM_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

//#ifndef CGAL_PLANAR_MAP_2_H
//#include <CGAL/Planar_map_2.h>
//#endif

#ifndef CGAL_INVERSE_INDEX_H
#include <CGAL/Inverse_index.h>
#endif

#include <iostream>

CGAL_BEGIN_NAMESPACE


template <class PM, class Drawer, class Window>
void draw_pm(const PM& pm,
             Drawer& drawer, 
             Window& w) 
{

  drawer.draw_faces(pm.faces_begin(), pm.faces_end());
  
  drawer.draw_halfedges(pm.halfedges_begin(), pm.halfedges_end());

  drawer.draw_vertices(pm.vertices_begin(), pm.vertices_end());
}

CGAL_END_NAMESPACE

#endif




