// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_UNIFORM_SIZING_FIELD_H
#define CGAL_MESH_3_UNIFORM_SIZING_FIELD_H

#include <CGAL/license/Mesh_3.h>


namespace CGAL {

namespace Mesh_3 {
  
template <typename Tr>
class Uniform_sizing_field
{
  typedef typename Tr::Geom_traits    Gt;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Gt::FT             FT;
  
public:
  // Vertices of mesh triangulation do not need to be updated 
  static const bool is_vertex_update_needed = false;
  
public:
  Uniform_sizing_field(const Tr&) {}
  void fill(const std::map<Weighted_point,FT>&) {}
  
  FT operator()(const Weighted_point&) const { return FT(1); }
  template <typename Handle>
  FT operator()(const Weighted_point&, const Handle&) const { return FT(1); }
};
  
  
} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_UNIFORM_SIZING_FIELD_H
