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
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/Alpha_shape_vertex_base_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef ALPHA_SHAPE_VERTEX_BASE_3_H
#define ALPHA_SHAPE_VERTEX_BASE_3_H

#include <utility>
#include <CGAL/Triangulation_vertex_base_3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class Gt>
class Alpha_shape_vertex_base_3 : public Triangulation_vertex_base_3<Gt>
{
public:

  typedef typename Gt::Coord_type Coord_type;
  typedef std::pair< Coord_type, Coord_type > Interval2;  
  typedef typename Gt::Point Point;

  //-------------------------------------------------------------------
private:

  Interval2 I;

  //-------------------------------------------------------------------
public:

  Alpha_shape_vertex_base_3()
    : Triangulation_vertex_base_3<Gt>()
    {}
  
  Alpha_shape_vertex_base_3(const Point& p)
    : Triangulation_vertex_base_3<Gt>(p)
    {}
  
  Alpha_shape_vertex_base_3(const Point& p, void* f)
    : Triangulation_vertex_base_3<Gt>(p, f) 
    {}

  //-------------------------------------------------------------------

  inline Interval2 get_range()
    {
      return I;
    }

  inline void set_range(Interval2 Inter)
    {
      I = Inter;
    }

};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //ALPHA_SHAPE_VERTEX_BASE_3_H
