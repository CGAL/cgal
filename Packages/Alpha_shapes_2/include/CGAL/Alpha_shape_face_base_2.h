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
// file          : include/CGAL/Alpha_shape_face_base_2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_ALPHA_SHAPE_FACE_BASE_2_H
#define CGAL_ALPHA_SHAPE_FACE_BASE_2_H

#include <CGAL/utility.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Df = Triangulation_face_base_2<Gt> >
class Alpha_shape_face_base_2 : public Df
{
  typedef Df                                         Base;
  typedef typename Df::Triangulation_data_structure  TDS;
public:
  typedef TDS                                Triangulation_data_structure;
  typedef typename TDS::Vertex_handle        Vertex_handle;
  typedef typename TDS::Face_handle          Face_handle;

  typedef typename Gt::Coord_type                    Coord_type;
  typedef Triple<Coord_type, Coord_type, Coord_type> Interval_3;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Base::template Rebind_TDS<TDS2>::Other    Fb2;
    typedef Alpha_shape_face_base_2<Gt,Fb2>         Other;
  };


private:
  Interval_3 vec_edge[3];
  Coord_type A;
 
public:
  Alpha_shape_face_base_2()  : Base()     {}
  
  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : Df(v0, v1, v2)     {}
  
  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
			  Face_handle n0, Face_handle n1, Face_handle n2)
    : Df(v0, v1, v2, n0, n1, n2)
    {}

  const Coord_type & get_alpha() const
    {
      return A;
    }

  void set_alpha(const Coord_type & AA)
    {
      A = AA;
    }

  const Interval_3 & get_ranges(int i) const
    {
      return vec_edge[i];
    }

  void set_ranges(int i, const Interval_3& Inter)
    {
      vec_edge[i]=Inter;
    }


};

CGAL_END_NAMESPACE

#endif //ALPHA_SHAPE_FACE_BASE_2_H
