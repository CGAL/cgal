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
// release       : $CGAL_Revision: CGAL-2.4-I-33 $
// release_date  : $CGAL_Date: 2001/12/04 $
//
// file          : include/CGAL/Topological_map_new.h
// package       : Planar_map (5.77)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
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




