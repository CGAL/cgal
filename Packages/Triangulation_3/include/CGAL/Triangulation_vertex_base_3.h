// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_base_vertex.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_Triangulation_BASE_VERTEX_H
#define CGAL_Triangulation_BASE_VERTEX_H

#include <CGAL/Triangulation_short_names.h>

template < class GT >
class CGAL_Triangulation_base_vertex 
{

public:
  typedef typename GT::Point Point;
  
  // CONSTRUCTORS
  
  inline 
  CGAL_Triangulation_base_vertex ()
    : _p(), _c(NULL)
  {}
  
  inline 
  CGAL_Triangulation_base_vertex(const Point & p)
    :  _p(p), _c(NULL)
  {}
    
  inline 
  CGAL_Triangulation_base_vertex(const Point & p, void* c)
    :  _p(p), _c(c)
  {}

  // ACCES 

  inline 
  Point point() const
  { return _p; }
    
  inline 
  void* cell() const
  { return _c; }

  // SETTING

  inline 
  void set_point(const Point & p)
  { _p = p; }
    
  inline 
  void set_cell(void* c)
  { _c = c; }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes 
  // to add their own purpose checking
  bool is_valid() const
  { return true; }


private:
  Point _p;
  void * _c;
  
};

#endif CGAL_TETRAHEDRALIZATION_BASE_VERTEX_H
