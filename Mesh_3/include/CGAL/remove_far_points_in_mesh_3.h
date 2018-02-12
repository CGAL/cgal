// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Clement Jamin
//
//******************************************************************************
// File Description : remove_far_points_in_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_REMOVE_FAR_POINTS_IN_MESH_3_H
#define CGAL_REMOVE_FAR_POINTS_IN_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {
  
namespace Mesh_3 {

/************************************************
// Class Mesher_3_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename C3T3, typename Concurrency_tag = typename C3T3::Concurrency_tag>
class Remove_far_points
{
#ifdef CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE

public:
  Remove_far_points(C3T3 &c3t3) : m_c3t3(c3t3) {}
  void remove_far_points() { m_c3t3.remove_far_points(); }
private:
  C3T3 &m_c3t3;

#else // !CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE

public:
  Remove_far_points(C3T3 &) {}
  void remove_far_points() {}

#endif
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename C3T3>
class Remove_far_points<C3T3, Parallel_tag>
{
public:
  Remove_far_points(C3T3 &c3t3) : m_c3t3(c3t3) {}
  
  void remove_far_points()
  {
    m_c3t3.remove_far_points();
  }

private:
  C3T3 &m_c3t3;
};
#endif // CGAL_LINKED_WITH_TBB

} // namespace Mesh_3

template <typename C3T3>
void
remove_far_points_in_mesh_3(C3T3& c3t3)
{
  typedef typename Mesh_3::Remove_far_points<C3T3> Remove_far_points;
  Remove_far_points cu(c3t3);
  cu.remove_far_points();
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REMOVE_FAR_POINTS_IN_MESH_3_H
