// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Julia
// extended class for regular_neighbor_coordinates computation 


#ifndef CGAL_REGULAR_NEIGHBOR_COORDINATES_TRAITS_2_H
#define CGAL_REGULAR_NEIGHBOR_COORDINATES_TRAITS_2_H

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>

CGAL_BEGIN_NAMESPACE 

template < class R, class W = typename  R::RT>
class Regular_neighbor_coordinates_traits_2
  : public Regular_triangulation_euclidean_traits_2<R,W>
{
public:
  typedef R                                     Rep;
  typedef typename R::FT                        FT;
  typedef typename Rep::Compute_area_2          Compute_area_2;
 
  Compute_area_2 compute_area_2_object() const
    {return Compute_area_2();}
};
 
CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_NEIGHBOR_COORDINATES_TRAITS_2_H
