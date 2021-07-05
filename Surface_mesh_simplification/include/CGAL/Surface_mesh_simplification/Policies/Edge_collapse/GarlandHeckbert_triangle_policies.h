// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_TRIANGLE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_TRIANGLE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_plane_quadrics.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_triangle_faces.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_optimizers.h>

#include <Eigen/Dense>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_simplification {

// use triangle quadrics for the faces, classical plane quadrics for the edges 
// and optimize with a check for singular matrices
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_triangle_policies : 
  public internal::GarlandHeckbert_triangle_faces<TriangleMesh, GeomTraits>,
  public internal::GarlandHeckbert_plane_edges<TriangleMesh, GeomTraits>,
  public internal::GarlandHeckbert_singular_optimizer<GeomTraits>,
  public internal::GarlandHeckbert_cost_base<
    typename boost::property_map<
      TriangleMesh, 
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_triangle_policies<TriangleMesh, GeomTraits>
  >,
  public internal::GarlandHeckbert_placement_base<
    typename boost::property_map<
      TriangleMesh, 
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_triangle_policies<TriangleMesh, GeomTraits>
  >
{
  public:

  typedef typename GeomTraits::FT FT;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;
  
  typedef CGAL::dynamic_vertex_property_t<Mat_4> Cost_property;
  typedef typename boost::property_map<TriangleMesh, Cost_property>::type Vertex_cost_map;
  
  typedef internal::GarlandHeckbert_placement_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_triangle_policies<TriangleMesh, GeomTraits>
    > Placement_base;

  typedef internal::GarlandHeckbert_cost_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_triangle_policies<TriangleMesh, GeomTraits>
    > Cost_base;
  
  using Cost_base::operator();
  using Placement_base::operator();
  
  GarlandHeckbert_triangle_policies(TriangleMesh& tmesh, FT dm = FT(100))
  {
    vcm_ = get(Cost_property(), tmesh);

    Cost_base::init_vcm(vcm_);
    Placement_base::init_vcm(vcm_);
  }

  typedef GarlandHeckbert_triangle_policies Get_cost;
  typedef GarlandHeckbert_triangle_policies Get_placement;
  
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }
  
  private:
  
  Vertex_cost_map vcm_;
};

} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
