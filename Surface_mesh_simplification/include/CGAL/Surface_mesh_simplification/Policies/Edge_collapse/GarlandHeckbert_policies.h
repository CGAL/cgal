// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <CGAL/tags.h>

#include <Eigen/Dense>

#include <boost/optional/optional.hpp>

#include <utility>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template <typename VCM_, typename FT_>
class GarlandHeckbert_cost
{
public:
  typedef VCM_                                                     Vertex_cost_map;
  typedef FT_                                                      FT;

  // Tells the main function of `Edge_collapse` that these policies must call "initialize"
  // and "update" functions. A bit awkward, but still better than abusing visitors.
  typedef CGAL::Tag_true                                           Update_tag;

  GarlandHeckbert_cost() { }
  GarlandHeckbert_cost(Vertex_cost_map vcm,
                       const FT discontinuity_multiplier = FT(100)) // abusing FT(double)
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
  FT m_discontinuity_multiplier;
};

template<typename VCM_>
class GarlandHeckbert_placement
{
public:
  typedef VCM_                                                                             Vertex_cost_map;

  GarlandHeckbert_placement() { }
  GarlandHeckbert_placement(Vertex_cost_map cost_matrices)
    : m_cost_matrices(cost_matrices)
  { }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile) const
  {
    typedef typename Profile::Triangle_mesh                                                Triangle_mesh;
    typedef typename Profile::Vertex_point_map                                             Vertex_point_map;
    typedef typename Profile::Geom_traits                                                  Geom_traits;
    typedef internal::GarlandHeckbert_core<Triangle_mesh, Vertex_point_map, Geom_traits>   GH_core;

    typedef typename GH_core::Matrix4x4                                                    Matrix4x4;
    typedef typename GH_core::Col4                                                         Col4;

    CGAL_precondition(get(m_cost_matrices, profile.v0()) != Matrix4x4());
    CGAL_precondition(get(m_cost_matrices, profile.v1()) != Matrix4x4());

    // the combined matrix has already been computed in the evaluation of the cost...
    const Matrix4x4 combinedMatrix = GH_core::combine_matrices(
                                       get(m_cost_matrices, profile.v0()),
                                       get(m_cost_matrices, profile.v1()));

    const Col4 p0 = GH_core::point_to_homogenous_column(profile.p0());
    const Col4 p1 = GH_core::point_to_homogenous_column(profile.p1());
    const Col4 opt = GH_core::optimal_point(combinedMatrix, p0, p1);

    boost::optional<typename Profile::Point> pt = typename Profile::Point(opt(0) / opt(3),
                                                                          opt(1) / opt(3),
                                                                          opt(2) / opt(3));

    return pt;
  }

private:
  Vertex_cost_map m_cost_matrices;
};

} // namespace internal

template <typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_policies
{
public:
  typedef TriangleMesh                                                          Triangle_mesh;
  typedef GeomTraits                                                            Geom_traits;
  typedef typename Geom_traits::FT                                              FT;

  typedef typename internal::GarlandHeckbert_matrix_type<GeomTraits>::type      Cost_matrix;
  typedef CGAL::dynamic_vertex_property_t<Cost_matrix>                          Cost_property;
  typedef typename boost::property_map<TriangleMesh, Cost_property>::type       Vertex_cost_map;

  typedef internal::GarlandHeckbert_cost<Vertex_cost_map, FT>                   Get_cost;
  typedef internal::GarlandHeckbert_placement<Vertex_cost_map>                  Get_placement;

  GarlandHeckbert_policies(TriangleMesh& tmesh,
                           const FT discontinuity_multiplier = FT(100))
  {
    vcm_ = get(Cost_property(), tmesh);
    get_cost_ = Get_cost(vcm_, discontinuity_multiplier);
    get_placement_ = Get_placement(vcm_);
  }

  Get_cost& get_cost() { return get_cost_; }
  const Get_cost& get_cost() const { return get_cost_; }
  Get_placement& get_placement() { return get_placement_; }
  const Get_placement& get_placement() const { return get_placement_; }

private:
  Vertex_cost_map vcm_;
  Get_cost get_cost_;
  Get_placement get_placement_;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
