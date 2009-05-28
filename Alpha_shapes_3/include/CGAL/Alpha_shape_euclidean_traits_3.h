// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

CGAL_BEGIN_NAMESPACE

template <class K>
class Alpha_shape_euclidean_traits_3 : public K
{
public:

 class  Compute_squared_radius_3 {
    typedef typename K::Point_3                  Point_3;
    typedef typename K::FT                       FT;
    typedef typename K::Compute_squared_radius_3   Compute_squared_radius_base;
 
  public:
    FT operator() (Point_3 p, 
		   Point_3 q , 
		   Point_3 r, 
		   Point_3 s) {
      return Compute_squared_radius_base()(p,q,r,s); }

    FT operator() (Point_3 p, 
		   Point_3 q , 
		   Point_3 r) {
      return Compute_squared_radius_base()(p,q,r); }

    FT operator() (Point_3 p, 
		   Point_3 q ) {
      return Compute_squared_radius_base()(p,q); }

    FT operator() (Point_3 p) {
      return FT(0);}
 };

//---------------------------------------------------------------------

  Compute_squared_radius_3 
  compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_3();
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
