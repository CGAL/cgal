// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : Triangulation/include/CGAL/Regular_triangulation_vertex_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/Triangulation_vertex_base_2.h>

CGAL_BEGIN_NAMESPACE

template < class GT >
class Regular_triangulation_vertex_base_2 : 
  public Triangulation_vertex_base_2<GT> {

public:
  typedef typename GT::Point_2 Point;

  Regular_triangulation_vertex_base_2 ()
    : Triangulation_vertex_base_2<GT>(), _hidden(false)
    {}
    
  Regular_triangulation_vertex_base_2(const Point & p, void * f = NULL)
    :  Triangulation_vertex_base_2<GT>(p, f), _hidden(false)
    {}

  void set_hidden(bool b) { _hidden = b; }
  bool is_hidden() { return _hidden ;}
 
private:
  bool _hidden;

};

CGAL_END_NAMESPACE

#endif //CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
