// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H
#define CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H

#include <CGAL/Polyhedron_items_3.h>

namespace CGAL {

template <class Refs, class TriangulatedMixedComplex3>
struct Skin_Surface_polyhedral_face : public CGAL::HalfedgeDS_face_base<Refs> {
  typedef typename TriangulatedMixedComplex3::Cell_handle Triang_Cell_handle;

  Triang_Cell_handle triang_ch;
};

template < class TriangulatedMixedComplex3 >
class Skin_surface_polyhedral_items_with_face_information_3
  : public Polyhedron_items_3 {
  
  template <class Refs, class Traits>
  struct Face_wrapper {
    typedef Skin_Surface_polyhedral_face<Refs, TriangulatedMixedComplex3> Face;
  };
};

} //namespace CGAL

#endif // CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H
