// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SELF_INTERSECTION_EXCEPTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SELF_INTERSECTION_EXCEPTION_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <stdexcept>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {

struct Self_intersection_exception :
  public std::runtime_error
{
  Self_intersection_exception()
    : std::runtime_error("Self-intersection detected in input mesh")
  {}
};

} } } // end of CGAL::Polygon_mesh_processing::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SELF_INTERSECTION_EXCEPTION_H
