// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Baskin Burak Senbaslar
//

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H

#include <CGAL/license/Surface_mesh_simplification.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template <class FT_>
class GarlandHeckbert_cost_stop_predicate
{
  typedef FT_                                           FT;

public:
  GarlandHeckbert_cost_stop_predicate(const FT gh_cost_threshold)
    :
      m_gh_cost_threshold(gh_cost_threshold)
  { }

  template <typename Profile>
  bool operator()(const FT& current_cost,
                  const Profile&,
                  std::size_t,
                  std::size_t) const
  {
    return current_cost >= m_gh_cost_threshold;
  }

private:
  const FT m_gh_cost_threshold;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
