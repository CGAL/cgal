// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_COMPLEX_2_IN_TRIANGULATION_3_H
#define CGAL_SURFACE_MESH_COMPLEX_2_IN_TRIANGULATION_3_H

#include <CGAL/Complex_2_in_triangulation_3.h>

namespace CGAL {

/**  Simple forward class for CGAL::Complex_2_in_triangulation_3<Tr>.
 *
 *   CGAL::Surface_mesher_complex_2_in_triangulation_3<Tr> is an alias
 *   for CGAL::Complex_2_in_triangulation_3<Tr>.
 */
template <class Tr>
class Surface_mesh_complex_2_in_triangulation_3 :
    public Complex_2_in_triangulation_3<Tr>
{
public:
  Surface_mesh_complex_2_in_triangulation_3 (Tr& tr) 
    : Complex_2_in_triangulation_3<Tr>(tr)
  {
  }
}; // end Surface_mesh_complex_2_in_triangulation_3

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_COMPLEX_2_IN_TRIANGULATION_3_H
