// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : remove_isolated_vertices_in_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H
#define CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {

template <typename C3T3>
void
remove_isolated_vertices_in_mesh_3(C3T3& c3t3)
{
  c3t3.remove_isolated_vertices();
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H
