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

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <iostream>
#include <utility>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_>
struct GarlandHeckbert_edge_collapse_visitor_base : Edge_collapse_visitor_base<TM_>
{
  typedef TM_                                                                     TM;
  typedef typename Edge_collapse_visitor_base<TM>::vertex_descriptor              vertex_descriptor;
  typedef typename Edge_collapse_visitor_base<TM>::halfedge_descriptor            halfedge_descriptor;

  typedef typename Edge_collapse_visitor_base<TM>::Kernel                         Kernel;
  typedef typename internal::GarlandHeckbertCore<TM>::FT                          FT;

  typedef typename Edge_collapse_visitor_base<TM>::Profile                        Profile;
  typedef typename internal::GarlandHeckbertCore<TM>::Matrix4x4                   Matrix4x4;
  typedef typename internal::GarlandHeckbertCore<TM>::garland_heckbert_state_type garland_heckbert_state_type;

  GarlandHeckbert_edge_collapse_visitor_base(garland_heckbert_state_type& aCostMatrices,
                                             FT aDiscontinuityMultiplier = 100.0)
    :
      mCostMatrices(aCostMatrices),
      mDiscontinuityMultiplier(aDiscontinuityMultiplier)
  { }

  void OnStarted(TM& aTM)
  {
    for(halfedge_descriptor& hd: halfedges(aTM))
    {
      vertex_descriptor target_vertex = target(hd, aTM);
      if(mCostMatrices.find(target_vertex) == mCostMatrices.end())
      {
        mCostMatrices.emplace(target_vertex,
          std::move(internal::GarlandHeckbertCore<TM>::fundamental_error_quadric(hd, aTM, mDiscontinuityMultiplier)));
      }
    }
  }

  void OnCollapsed(const Profile& aProfile,
                   const vertex_descriptor& aVD)
  {
    Matrix4x4 combinedMatrix =
        std::move(internal::GarlandHeckbertCore<TM>::combine_matrices(
                    mCostMatrices.at(aProfile.v0()),
                    mCostMatrices.at(aProfile.v1())));

    mCostMatrices.erase(aProfile.v0());
    mCostMatrices.erase(aProfile.v1());

    mCostMatrices.emplace(aVD, std::move(combinedMatrix));
  }

private:
  garland_heckbert_state_type& mCostMatrices;
  FT mDiscontinuityMultiplier;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_GARLANDHECKBERT_EDGE_COLLAPSE_VISITOR_BASE_H
