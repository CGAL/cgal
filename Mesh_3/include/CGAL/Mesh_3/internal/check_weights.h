// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H
#define CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/enum.h>
#include <CGAL/tags.h>
#include <CGAL/type_traits.h>
#include <CGAL/STL_Extension/internal/Has_features.h>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

template<typename Triangulation, typename MeshDomain>
bool has_non_protecting_weights(const Triangulation& tr, const MeshDomain&, Tag_true)
{
  constexpr bool with_features = ::CGAL::internal::Has_features<MeshDomain>::value;

  typedef typename Triangulation::FT                FT;

  auto cwsr = tr.geom_traits().compare_weighted_squared_radius_3_object();

  for (typename Triangulation::Finite_vertices_iterator
        vv = tr.finite_vertices_begin();
        vv != tr.finite_vertices_end();
        ++vv)
  {
    const auto& vv_wp = tr.point(vv);
    if (cwsr(vv_wp, FT(0)) != CGAL::EQUAL)
    {
      if constexpr (with_features)
      {
        if (vv->in_dimension() > 1)
          return true;
      }
      else
        return true;
    }
  }
  return false;
}

template <typename Triangulation, typename MeshDomain>
bool has_non_protecting_weights(const Triangulation&, const MeshDomain&, Tag_false)
{
  return false;
}

template<typename Triangulation, typename MeshDomain>
bool has_non_protecting_weights(const Triangulation& tr, const MeshDomain& md)
{
  return has_non_protecting_weights(tr, md, Boolean_tag<is_regular_triangulation_v<Triangulation> >());
}

} // end namespace internal
} // end namespace Mesh_3
} // end namespace CGAL

#endif //CGAL_INTERNAL_MESH_3_CHECK_WEIGHTS_H

