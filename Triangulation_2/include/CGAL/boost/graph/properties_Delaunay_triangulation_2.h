// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PROPERTIES_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PROPERTIES_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Delaunay_triangulation_2.h>

#define CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS typename GT, typename TDS
#define CGAL_2D_TRIANGULATION CGAL::Delaunay_triangulation_2<GT, TDS>

#include <CGAL/boost/graph/internal/properties_2D_triangulation.h>

#endif /* CGAL_PROPERTIES_DELAUNAY_TRIANGULATION_2_H */
