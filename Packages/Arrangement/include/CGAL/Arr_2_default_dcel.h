// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : ?/include/CGAL/Arr_2_default_dcel.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_2_DEFAULT_DCEL_H
#define CGAL_ARR_2_DEFAULT_DCEL_H

#ifndef CGAL_PM_DEFAULT_DCEL_H
#include <CGAL/Pm_default_dcel.h>
#endif

#ifndef ARR_2_BASES_H
#include <CGAL/Arr_2_bases.h>
#endif

CGAL_BEGIN_NAMESPACE


///////////////////////////////////////////////////////////////
//               DEFAULT DCEL
///////////////////////////////////////////////////////////////

template <class Traits>
class Arr_2_default_dcel
  : public Pm_dcel<
Arr_2_vertex_base<typename Traits::Point>,
Arr_2_halfedge_base<Arr_base_node<typename Traits::X_curve> >,
Arr_2_face_base
> 
{
public:  // CREATION
  
  Arr_2_default_dcel() {}
  
};



CGAL_END_NAMESPACE

#endif










