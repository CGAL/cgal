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
// file          : include/CGAL/Alpha_shape_vertex_base_2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_ALPHA_SHAPE_VERTEX_BASE_2_H
#define CGAL_ALPHA_SHAPE_VERTEX_BASE_2_H

#include <utility>
#include <CGAL/Triangulation_vertex_base_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class Gt, class Vb = Triangulation_vertex_base_2<Gt> >
class Alpha_shape_vertex_base_2 : public Vb
{
  typedef Vb                                         Base;
  typedef typename Vb::Triangulation_data_structure  TDS;
public:
  typedef TDS                             Triangulation_data_structure;
  typedef typename TDS::Vertex_handle     Vertex_handle;
  typedef typename TDS::Face_handle       Face_handle;

  typedef typename Gt::Coord_type                Coord_type;
  typedef std::pair< Coord_type, Coord_type >    Interval2;
  typedef typename Base::Point                   Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Base::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef Alpha_shape_vertex_base_2 <Gt,Vb2>         Other;
  };
private:
  Interval2 I;

public:
  Alpha_shape_vertex_base_2()
    : Base() 
    {}
  
  Alpha_shape_vertex_base_2(const Point & p)
    : Base(p) 
    {}
  
  Alpha_shape_vertex_base_2(const Point & p, Face_handle f)
    : Base(p, f) 
    {}


public:

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

#endif //ALPHA_SHAPE_VERTEX_BASE_2_H
