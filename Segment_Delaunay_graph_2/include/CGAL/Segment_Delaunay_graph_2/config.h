// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONFIG_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONFIG_H 1

#include <CGAL/basic.h>
#include <CGAL/version.h>

#if (CGAL_VERSION_NR < 1030210000)
#error "Incompatible CGAL version"
#elif ((CGAL_VERSION_NR > 1030210000) && (CGAL_VERSION_NR < 1030300000))
#define CGAL_CFG_NO_OPERATOR_TIMES_FOR_SIGN 1
#endif

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONFIG_H
