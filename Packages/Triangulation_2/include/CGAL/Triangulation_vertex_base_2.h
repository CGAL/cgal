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
// source        : $Source$
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

template < class GT >
class CGAL_Triangulation_vertex_base_2 {

public:
    typedef typename GT::Point Point;


  CGAL_Triangulation_vertex_base_2 ()
        : _p(),_f(NULL)
    {}
    
    CGAL_Triangulation_vertex_base_2(const Point & p)
        :  _p(p), _f(NULL)
    {}
    
    CGAL_Triangulation_vertex_base_2(const Point & p, void* f)
        :  _p(p), _f(f)
    {}

    inline void set_point(const Point & p)
    {
        _p = p;
    }
    
    inline void set_face(void* f)
    {
        _f = f;
    }

    inline 
    const Point&  point() const
    {
        return _p;
    }
    
    
    inline void* face() const
    {
        return _f;
    }
    
  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid() const
  {return true;}

private:
        Point _p;
        void * _f;

};

#endif CGAL_TRIANGULATION_VERTEX_BASE_2_H
