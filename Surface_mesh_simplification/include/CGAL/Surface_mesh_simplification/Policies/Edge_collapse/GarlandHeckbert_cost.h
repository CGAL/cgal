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
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <CGAL/tags.h>

#include <boost/optional/optional.hpp>

#include <utility>

namespace CGAL {
namespace Surface_mesh_simplification {

template <typename GT_>
struct GarlandHeckbert_cost_matrix
{
  typedef typename Eigen::Matrix<typename GT_::FT, 4, 4>           type;
};

template <typename VCM_>
class GarlandHeckbert_cost
{
public:
  typedef VCM_                                                     Vertex_cost_map;
  typedef typename boost::property_traits<VCM_>::value_type        Cost_matrix; // @tmp name

  // Tells the edge collapse main function that we need to call "initialize"
  // and "update" functions. A bit awkward, but still better than abusing visitors.
  typedef CGAL::Tag_true                                           Update_tag;

  GarlandHeckbert_cost(Vertex_cost_map& vcm,
                       const double discontinuity_multiplier = 100.) // abusing FT(double)
    :
      m_cost_matrices(vcm),
      m_discontinuity_multiplier(discontinuity_multiplier)
  { }

  template <typename TM_, typename VPM_, typename GT_>
  void initialize(const TM_& tmesh, const VPM_& vpm, const GT_& gt) const
  {
    typedef internal::GarlandHeckbert_core<TM_, VPM_, GT_>                                 GH_core;

    GH_core::fundamental_error_quadrics(m_cost_matrices, tmesh, vpm, gt, m_discontinuity_multiplier);
  }

  template <typename Profile>
  boost::optional<typename Profile::FT>
  operator()(const Profile& profile,
             const boost::optional<typename Profile::Point>& placement) const
  {
    typedef typename Profile::Triangle_mesh                                                Triangle_mesh;
    typedef typename Profile::Vertex_point_map                                             Vertex_point_map;
    typedef typename Profile::Geom_traits                                                  Geom_traits;
    typedef internal::GarlandHeckbert_core<Triangle_mesh, Vertex_point_map, Geom_traits>   GH_core;

    typedef typename GH_core::Matrix4x4                                                    Matrix4x4;
    typedef typename GH_core::Col4                                                         Col4;
    typedef boost::optional<typename Profile::FT>                                          Optional_FT;

    if(!placement)
      return boost::optional<typename Profile::FT>();

    CGAL_precondition(get(m_cost_matrices, profile.v0()) != Matrix4x4());
    CGAL_precondition(get(m_cost_matrices, profile.v1()) != Matrix4x4());

    const Matrix4x4 combined_matrix = GH_core::combine_matrices(get(m_cost_matrices, profile.v0()),
                                                                get(m_cost_matrices, profile.v1()));
    const Col4 pt = GH_core::point_to_homogenous_column(*placement);
    const Optional_FT cost = (pt.transpose() * combined_matrix * pt)(0, 0);

    return cost;
  }

  template <typename Profile, typename vertex_descriptor>
  void update_after_collapse(const Profile& profile,
                             const vertex_descriptor new_v) const
  {
    typedef typename Profile::Triangle_mesh                                                Triangle_mesh;
    typedef typename Profile::Vertex_point_map                                             Vertex_point_map;
    typedef typename Profile::Geom_traits                                                  Geom_traits;
    typedef internal::GarlandHeckbert_core<Triangle_mesh, Vertex_point_map, Geom_traits>   GH_core;

    put(m_cost_matrices, new_v,
        GH_core::combine_matrices(get(m_cost_matrices, profile.v0()),
                                  get(m_cost_matrices, profile.v1())));
  }

private:
  Vertex_cost_map m_cost_matrices;
  const double m_discontinuity_multiplier;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_H
