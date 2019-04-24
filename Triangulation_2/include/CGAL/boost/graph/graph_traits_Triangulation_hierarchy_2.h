// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
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
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H

#include <CGAL/Triangulation_hierarchy_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Delaunay_triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Constrained_triangulation_2 object.

#define CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS typename Tr
#define CGAL_2D_TRIANGULATION CGAL::Triangulation_hierarchy_2<Tr>
#define CGAL_2D_TRIANGULATION_TEMPLATES Tr

#include <CGAL/boost/graph/graph_traits_2D_triangulation.h>

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H
