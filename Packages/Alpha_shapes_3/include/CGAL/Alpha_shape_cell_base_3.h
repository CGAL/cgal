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
// file          : include/CGAL/Alpha_shape_cell_base_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef ALPHA_SHAPE_CELL_BASE_3_H
#define ALPHA_SHAPE_CELL_BASE_3_H

#include <vector>
#include <CGAL/triple.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Gt, class Df >
class Alpha_shape_cell_base_3 : public Df
{
public:

  typedef typename Gt::Rep Rp;               // A corriger avec Monique
  typedef typename Rp::FT Coord_type;        // A corriger avec Monique
  typedef triple<Coord_type, Coord_type , Coord_type> Interval3;

  //-------------------------------------------------------------------
private:

  Interval3 vec_facet[4];
  Interval3 vec_edge[4][4];
  Coord_type A;

  //-------------------------------------------------------------------
public:
  
  Alpha_shape_cell_base_3() 
    : Df()
    {}
  
  Alpha_shape_cell_base_3(void* v0, void* v1, void* v2, void* v3)
    : Df( v0, v1, v2, v3)
    {}
  
  Alpha_shape_cell_base_3(void* v0, void* v1, void* v2, void* v3,
			  void* n0, void* n1, void* n2, void* n3)
    : Df(v0,  v1,  v2, v3,
	 n0,  n1,  n2, n3)
    {}

  //-------------------------------------------------------------------

  inline Coord_type get_alpha()
    {
      return A;
    }
  
  inline void set_alpha(Coord_type AA)
    {
      A = AA;
    }

  inline Interval3 get_facet_ranges(const int& i)
    {
      return vec_facet[i];
    }

  inline void set_facet_ranges(const int& i, const Interval3& Inter)
    {
      vec_facet[i]=Inter;
    }
  
  inline Interval3 get_edge_ranges(const int& i, const int& j)
    {
      return vec_edge[i][j];
    }

  inline void set_edge_ranges(const int& i, const int& j, 
			      const Interval3& Inter)
    {
      vec_edge[i][j]=Inter;
    }
  
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //ALPHA_SHAPE_CELL_BASE_3_H
