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
// file          : Triangulation/include/CGAL/Triangulation_vertex_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class GT >
class Triangulation_vertex_base_2 {

public:
  typedef typename GT::Point_2 Point;

  Triangulation_vertex_base_2 ()
    : _p(Point()), _f(NULL)
    {}
    
  Triangulation_vertex_base_2(const Point & p, void * f = NULL)
    :  _p(p), _f(f)
    {}

  void set_point(const Point & p) { _p = p; }
  void set_face(void* f) { _f = f ;}
 
  const Point&  point() const { return _p; }
  void* face() const { return _f;}
 
  // the non const version of point() is undocument
  // but needed to make the point iterator works
  // using Lutz projection scheme
  Point&        point() { return _p; }

    
  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
    {return true;}

private:
        Point _p;
        void * _f;

};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_VERTEX_BASE_2_H
