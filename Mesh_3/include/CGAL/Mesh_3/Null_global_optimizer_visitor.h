// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_NULL_GLOBAL_OPTIMIZER_VISITOR_H
#define CGAL_MESH_3_NULL_GLOBAL_OPTIMIZER_VISITOR_H

#include <CGAL/license/Mesh_3.h>


namespace CGAL {
namespace Mesh_3 {

template < typename C3T3 >
class Null_global_optimizer_visitor
{
  typedef typename C3T3::Triangulation    Tr;
  typedef typename Tr::Geom_traits::FT    FT;

public:
  void after_compute_moves() {}
  void after_move_points() {}
  void after_rebuild_restricted_delaunay() {}
  void end_of_iteration(int) {}
};

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_NULL_GLOBAL_OPTIMIZER_VISITOR_H
