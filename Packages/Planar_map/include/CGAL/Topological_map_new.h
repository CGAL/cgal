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
#ifndef  CGAL_TOPOLOGICAL_MAP_H
#define  CGAL_TOPOLOGICAL_MAP_H

#include <CGAL/Topological_map_items.h>
#include <CGAL/HalfedgeDS_default.h>

CGAL_BEGIN_NAMESPACE

template < class Traits,
           class TopologicalMapItems = Topological_map_items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class T_HDS = HalfedgeDS_default, 
           class Alloc = CGAL_ALLOCATOR(int)>
class Topological_map
{
  typedef Topological_map < Traits, TopologicalMapItems, T_HDS, Alloc> Self;
  typedef TopologicalMapItems                                          Items;

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  typedef T_HDS< Traits, Items, Alloc>  HDS;
#else
  typedef typename T_HDS::template HDS< Traits, Items, Alloc>  HDS;
#endif

};

CGAL_END_NAMESPACE

#endif // CGAL_TOPOLOGICAL_MAP_H




