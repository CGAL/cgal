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
// file          : include/CGAL/Triangulation_vertex_base_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_short_names_3.h>

template < class GT >
class CGAL_Triangulation_vertex_base_3
{

public:
  typedef typename GT::Point Point;
  
  // CONSTRUCTORS
  
  inline 
  CGAL_Triangulation_vertex_base_3()
    : _p(), _c(NULL)
  {}
  
  inline 
  CGAL_Triangulation_vertex_base_3(const Point & p)
    :  _p(p), _c(NULL)
  {}
    
  inline 
  CGAL_Triangulation_vertex_base_3(const Point & p, void* c)
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
  bool is_valid(bool verbose, int level) const
  { 
    //return true; 
    return( cell() != NULL );
  }


private:
  Point _p;
  void * _c;
  
};

#endif CGAL_TRIANGULATION_VERTEX_BASE_3_H
