// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PROPERTIES_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_PROPERTIES_CONSTRAINED_TRIANGULATION_2_H

#include <CGAL/Constrained_triangulation_2.h>

#define CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS typename GT, typename TDS, typename Itag
#define CGAL_2D_TRIANGULATION CGAL::Constrained_triangulation_2<GT, TDS, Itag>

#include <CGAL/boost/graph/properties_2D_triangulation.h>

#endif /* CGAL_PROPERTIES_CONSTRAINED_TRIANGULATION_2_H */
