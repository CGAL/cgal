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

template < class Gt, class Df >
class Alpha_shape_face_base_2 : public Df
{
public:
  typedef typename Gt::Coord_type Coord_type;
  typedef Triple<Coord_type, Coord_type, Coord_type> Interval_3;

  Alpha_shape_face_base_2() 
    : Df() 
    {}
  
  Alpha_shape_face_base_2(void* v0, void* v1, void* v2)
    : Df(v0, v1, v2)
    {}
  
  Alpha_shape_face_base_2(void* v0, void* v1, void* v2,
			  void* n0, void* n1, void* n2)
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

private:
  Interval_3 vec_edge[3];
  Coord_type A;
};

CGAL_END_NAMESPACE

#endif //ALPHA_SHAPE_FACE_BASE_2_H
