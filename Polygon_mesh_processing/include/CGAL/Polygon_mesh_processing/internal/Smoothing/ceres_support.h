// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

// Ignore warnings on Windows of type
// " class XXXX needs to have dll-interface to be used by clients of struct YYYY" in an interaction
// between Ceres and Eigen, which doesn't make much sense as Eigen is header-only.

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4251)
#endif

#include "ceres/ceres.h"

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H
