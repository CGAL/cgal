// ============================================================================
//
// Copyright (c) 1998, 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_hierarchy_vertex_base_3.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Triangulation3
// author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb>
class Triangulation_hierarchy_vertex_base_3
  : public Vbb
{
public:
  typedef Vbb                      V_Base;
  typedef typename V_Base::Point   Point;

  Triangulation_hierarchy_vertex_base_3()
    : V_Base(), _up(0), _down(0) {}

  Triangulation_hierarchy_vertex_base_3(const Point & p, void* f)
    : V_Base(p,f), _up(0), _down(0) {}

  Triangulation_hierarchy_vertex_base_3(const Point & p)
    : V_Base(p), _up(0), _down(0) {}

  void* up() const {return _up;}
  void* down() const {return _down;}
  void set_up(void *u) {_up=u;}
  void set_down(void *d) {if (this) _down=d;}

private:
  void* _up;    // same vertex one level above
  void* _down;  // same vertex one level below
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H
