// Copyright (c) 1999  Tel-Aviv University (Israel).
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
// Author(s)     : Iddo Hanniel
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
  Arr_2_halfedge_base<Arr_base_node<typename Traits::Curve, 
    typename Traits::X_curve> >, Arr_2_face_base > 
{
public:  // CREATION
  
  Arr_2_default_dcel() {}
  
};



CGAL_END_NAMESPACE

#endif










