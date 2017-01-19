// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Naama mayer         <naamamay@post.tau.ac.il>


#ifndef CGAL_ARR_POLYHEDRAL_SGM_TRANSFORMATION_H
#define CGAL_ARR_POLYHEDRAL_SGM_TRANSFORMATION_H

#include <CGAL/license/Arrangement_on_surface_2.h>


namespace CGAL {

/* This function rotates the face when the arrangement is Arr_polyhedral_sgm */

template <class Arrangement, class Transformation_3>
class Arr_polyhedral_sgm_transformation
{
public:

  typedef typename Arrangement::Face_handle   Face_handle;

  void rotate_face(Face_handle f, const Transformation_3 & aff)
  {			
    //Transform all the vertices of the original polyhedron.
    f->set_point(aff.transform(f->point()));
  }
};

} //namespace CGAL

#endif
