// Copyright (C) 2018  GeometryFactory Sarl
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//

#ifndef CGAL_POLYHEDRON_3_FWD_H
#define CGAL_POLYHEDRON_3_FWD_H

#include <CGAL/memory.h>

/// \file Polyhedron_3_fwd.h
/// Forward declarations of the Polyhedron_3 package.

#ifndef DOXYGEN_RUNNING
namespace CGAL {

  class Polyhedron_items_3;

  template < class T, class I, class A>
  class HalfedgeDS_default;

  template < class PolyhedronTraits_3,
             class PolyhedronItems_3 = Polyhedron_items_3,
             template < class T, class I, class A>
             class T_HDS = HalfedgeDS_default, 
             class Alloc = CGAL_ALLOCATOR(int)
             >
  class Polyhedron_3;
  
} // CGAL
#endif

#endif /* CGAL_POLYHEDRON_3_FWD_H */


