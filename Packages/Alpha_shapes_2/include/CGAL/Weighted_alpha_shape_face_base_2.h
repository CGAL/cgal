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
// file          : include/CGAL/Weighted_alpha_shape_face_base_2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_FACE_BASE_2_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_FACE_BASE_2_H

#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/triple.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Gt >
class Weighted_alpha_shape_face_base_2 : public
Regular_triangulation_face_base_2<Gt> 
{
  
private:
  typedef typename Gt::Coord_type Coord_type;
  typedef triple<Coord_type, Coord_type , Coord_type> Interval3;
  Interval3 vec_edge[3];
  Coord_type A;

public:
  inline Coord_type get_Alpha() 
    {
      return A;
    }
  
  inline void Alpha(Coord_type AA) 
    {
      A = AA;
    }


  inline Interval3 get_Range(int i) 
    {
      return vec_edge[i];
    }

  inline void Range(int i, Interval3 Inter) 
    {
      vec_edge[i]=Inter;
    }
  
  
  Weighted_alpha_shape_face_base_2() 
    : Regular_triangulation_face_base_2<Gt>() 
    {}

  
  Weighted_alpha_shape_face_base_2(void* v0, void* v1, void* v2)
    : Regular_triangulation_face_base_2<Gt>( v0, v1, v2) 
    {}

  
  
  Weighted_alpha_shape_face_base_2(void* v0, void* v1, void* v2,
				   void* n0, void* n1, void* n2)
    :  Regular_triangulation_face_base_2<Gt>(v0,  v1,  v2,
					     n0,  n1,  n2) 
    {}
  
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_WEIGHTED_ALPHA_SHAPE_FACE_BASE_2_H
