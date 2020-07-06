// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CURVES_ON_SURFACE_TOPOLOGY_H
#define CGAL_CURVES_ON_SURFACE_TOPOLOGY_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Surface_mesh_topology/internal/Minimal_quadrangulation.h>
#include <CGAL/Surface_mesh_topology/internal/Shortest_noncontractible_cycle.h>
#include <CGAL/Surface_mesh_topology/internal/Facewidth.h>
#include <CGAL/Face_graph_wrapper.h>
#include <memory>

namespace CGAL {
namespace Surface_mesh_topology {

template<typename Mesh>
class Curves_on_surface_topology
{
public:
  // Types for minimal quadrangulation
  using Minimal_quadrangulation=internal::Minimal_quadrangulation<Mesh>;
  using Reduced_map            =typename Minimal_quadrangulation::Local_map;

  // Types for shortest noncontractible cycle
  using Shortest_noncontractible_cycle=typename internal::Shortest_noncontractible_cycle<Mesh, true>;
  using Facewidth                     =typename internal::Facewidth<Mesh>;
  using Dart_const_handle             =typename Shortest_noncontractible_cycle::Original_dart_const_handle;
  using halfedge_descriptor           =Dart_const_handle ; // To be compatible with BGL

  // Constructor
  Curves_on_surface_topology(const Mesh& amesh, bool /* display_time */=false) :
    m_original_mesh(amesh),
    m_minimal_quadrangulation(nullptr),
    m_shortest_noncontractible_cycle(nullptr),
    m_facewidth(nullptr)
  {}

//================================================================================
// Homotopy test

  /// @return true iff the minimal quadrangulation is computed.
  bool is_minimal_quadrangulation_computed() const
  { return m_minimal_quadrangulation!=nullptr; }

  /// Computes the minimal quadrangulation if it is not yet computed.
  void compute_minimal_quadrangulation(bool display_time=true) const
  {
    if (m_minimal_quadrangulation==nullptr)
    {
      m_minimal_quadrangulation=std::make_unique<Minimal_quadrangulation>
        (m_original_mesh, display_time);
    }
  }

  /// Return the reduced map computed in the minimal quadrangulation.
  /// @pre is_minimal_quadrangulation_computed()
  const Reduced_map& get_minimal_quadrangulation() const
  {
    CGAL_assertion(is_minimal_quadrangulation_computed());
    return m_minimal_quadrangulation->get_local_map();
  }

  /// @return true iff 'p' is contractible.
  bool is_contractible(const Path_on_surface<Mesh>& p,
                       bool display_time=false) const
  {
    compute_minimal_quadrangulation(display_time);
    return m_minimal_quadrangulation->is_contractible(p, display_time);
  }

  /// @return true iff 'p1' and 'p2' are freely homotopic.
  bool are_freely_homotopic(const Path_on_surface<Mesh>& p1,
                            const Path_on_surface<Mesh>& p2,
                            bool display_time=false) const
  {
    compute_minimal_quadrangulation(display_time);
    return m_minimal_quadrangulation->are_freely_homotopic(p1, p2,
                                                           display_time);
  }

  /// @return true iff 'p1' and 'p2' are base point freely homotopic.
  bool are_base_point_homotopic(const Path_on_surface<Mesh>& p1,
                                const Path_on_surface<Mesh>& p2,
                                bool display_time=false) const
  {
    compute_minimal_quadrangulation(display_time);
    return m_minimal_quadrangulation->are_base_point_homotopic(p1, p2,
                                                               display_time);
  }

//================================================================================
// Shortest noncontractible cycle; edge and face width

  bool is_shortest_noncontractible_cycle_representation_computed() const
  { return m_shortest_noncontractible_cycle!=nullptr; }

  void compute_shortest_non_contractible_cycle_representation(bool display_time=false) const
  {
    if (m_shortest_noncontractible_cycle==nullptr)
    {
      m_shortest_noncontractible_cycle=std::make_unique
        <Shortest_noncontractible_cycle>(m_original_mesh, display_time);
    }
  }

  template <class WeightFunctor>
  Path_on_surface<Mesh> compute_shortest_non_contractible_cycle_with_base_point
  (Dart_const_handle dh, const WeightFunctor& wf, bool display_time=false) const
  {
    compute_shortest_non_contractible_cycle_representation(display_time);
    return m_shortest_noncontractible_cycle->compute_cycle(dh, NULL, wf, display_time);
  }

  Path_on_surface<Mesh> compute_shortest_non_contractible_cycle_with_base_point
  (Dart_const_handle dh, bool display_time=false) const
  {
    compute_shortest_non_contractible_cycle_representation(display_time);
    return m_shortest_noncontractible_cycle->compute_cycle(dh, display_time);
  }

  Path_on_surface<Mesh> compute_shortest_non_contractible_cycle
  (bool display_time=false) const
  {
    compute_shortest_non_contractible_cycle_representation(display_time);
    return m_shortest_noncontractible_cycle->
      compute_shortest_non_contractible_cycle(nullptr, display_time);
  }

  template <class WeightFunctor>
  Path_on_surface<Mesh> compute_shortest_non_contractible_cycle
  (const WeightFunctor& wf, bool display_time=false) const
  {
    compute_shortest_non_contractible_cycle_representation(display_time);
    return m_shortest_noncontractible_cycle->
      compute_shortest_non_contractible_cycle(nullptr, wf, display_time);
  }

  Path_on_surface<Mesh> compute_edge_width(bool display_time=false) const
  {
    compute_shortest_non_contractible_cycle_representation(display_time);
    return m_shortest_noncontractible_cycle->compute_edge_width(display_time);
  }

  bool is_facewidth_representation_computed() const
  { return m_facewidth!=nullptr; }

  void compute_face_width_representation(bool display_time=false) const
  {
    if (m_facewidth==nullptr)
    { m_facewidth=std::make_unique<Facewidth>(m_original_mesh, display_time); }
  }

  std::vector<Dart_const_handle> compute_face_width(bool display_time=false) const
  {
    compute_face_width_representation(display_time);
    return m_facewidth->compute_face_width(display_time);
  }

protected:
  const Mesh&                                             m_original_mesh;
  mutable std::unique_ptr<Minimal_quadrangulation>        m_minimal_quadrangulation;
  mutable std::unique_ptr<Shortest_noncontractible_cycle> m_shortest_noncontractible_cycle;
  mutable std::unique_ptr<Facewidth>                      m_facewidth;
};

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_CURVES_ON_SURFACE_TOPOLOGY_H //
// EOF //
