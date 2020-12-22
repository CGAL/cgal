// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H

#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/boost/graph/properties_Triangulation_hierarchy_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Delaunay_triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Constrained_triangulation_2 object.

#define CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS typename Tr
#define CGAL_2D_TRIANGULATION CGAL::Triangulation_hierarchy_2<Tr>
#define CGAL_2D_TRIANGULATION_TEMPLATES Tr

#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation.h>

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_HIERARCHY_2_H
