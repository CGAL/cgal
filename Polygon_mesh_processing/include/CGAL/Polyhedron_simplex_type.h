// Copyright (c) 2014 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_POLYHEDRON_SIMPLEX_TYPE_H
#define CGAL_POLYHEDRON_SIMPLEX_TYPE_H

namespace CGAL{

enum Polyhedron_simplex_type {
    POLYHEDRON_VERTEX,
    POLYHEDRON_EDGE,
    POLYHEDRON_FACET,
    POLYHEDRON_NONE};
}

#endif // CGAL_POLYHEDRON_SIMPLEX_TYPE_H
