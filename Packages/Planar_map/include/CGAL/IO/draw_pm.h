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
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>
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




