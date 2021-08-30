// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/Lindstrom_Turk_core.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_>
class LindstromTurk_placement
{
public:
  typedef TM_                                                             TM;
  typedef internal::LindstromTurk_params                                  LindstromTurk_params;

  LindstromTurk_placement(const LindstromTurk_params& LT_params = LindstromTurk_params())
    : m_LT_params(LT_params)
  {}

  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& profile) const
  {
    return internal::LindstromTurkCore<TM,Profile>(m_LT_params, profile).compute_placement();
  }

private:
  LindstromTurk_params m_LT_params;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PLACEMENT_H
