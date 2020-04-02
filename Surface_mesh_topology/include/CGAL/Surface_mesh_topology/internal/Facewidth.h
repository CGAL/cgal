// Copyright (c) 2020 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Thien Hoang <thienvhoang99@gmail.com>
//
#ifndef CGAL_FACEWIDTH_H
#define CGAL_FACEWIDTH_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Path_on_surface.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_topology/internal/Generic_map_selector.h>
#include <CGAL/Surface_mesh_topology/internal/Shortest_noncontractible_cycle.h>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

template <class Mesh_>
class Facewidth
{
public:
  using Self=Facewidth<Mesh_>;
  using Mesh=Mesh_;

  using Original_map_wrapper=internal::Generic_map_selector<Mesh>;
  using Original_dart_const_handle=typename Original_map_wrapper::Dart_const_handle_original;

  using Local_map        =typename Original_map_wrapper::Generic_map;
  using Dart_handle      =typename Local_map::Dart_handle;
  using Dart_const_handle=typename Local_map::Dart_const_handle;
  using size_type        = typename Local_map::size_type;

  using Path             = CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;
  using SNC              = Shortest_noncontractible_cycle<Local_map>;

  Facewidth(const Mesh& amesh, bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    Local_map m_gmap; // TODO REMOVE
    typename Original_map_wrapper::Origin_to_copy_map origin_to_radial;
    Original_map_wrapper::copy(m_radial_map, const_cast<Mesh&>(amesh),
                               origin_to_radial, m_copy_to_origin, Local_map::NB_MARKS);
    m_copy_to_origin.clear();
    Original_map_wrapper::copy(m_gmap, const_cast<Mesh&>(amesh),
                               m_origin_to_copy, m_copy_to_origin, Local_map::NB_MARKS); // TODO BETTER
    // Initialize 0-attributes for m_gmap
    int counter = 0;
    for (auto it=m_gmap.darts().begin(), itend=m_gmap.darts().end(); it!=itend; ++it)
    {
      if (m_gmap.template attribute<0>(it)==NULL)
      {
        m_gmap.template set_attribute<0>(it, m_gmap.template create_attribute<0>(counter++));
        m_vertex_list.push_back(it);
      }
    }

    // m_face_list contains dart handles of m_gmap
    for (auto it=m_gmap.template one_dart_per_cell<2>().begin(),
         itend=m_gmap.template one_dart_per_cell<2>().end(); it!=itend; ++it)
    { m_face_list.push_back(it); }

    // Create edge_list
    std::vector<Dart_handle> edge_list;
    for (auto it=m_radial_map.template one_dart_per_cell<1>().begin(),
         itend=m_radial_map.template one_dart_per_cell<1>().end(); it!=itend; ++it)
    { edge_list.push_back(it); }

    // m_radial_map hasn't been changed so far

    // face_list contains dart handles of m_radial_map
    std::vector<Dart_handle> face_list;
    for (auto it=m_radial_map.template one_dart_per_cell<2>().begin(),
         itend=m_radial_map.template one_dart_per_cell<2>().end(); it!=itend; ++it)
    { face_list.push_back(it); }

    // Adding "centroids"
    std::vector<Dart_handle> centroids;
    for (auto it : face_list) // face_list contains dart handles of m_radial_map
    {
      auto new_vertex=m_radial_map.insert_cell_0_in_cell_2(it);
      centroids.push_back(new_vertex);
    }

    // Initialize 1-attributes of m_radial_map
    for (auto it=m_radial_map.darts().begin(), itend=m_radial_map.darts().end(); it!=itend; ++it)
    {
      if (m_radial_map.template attribute<1>(it)==NULL)
      { m_radial_map.template set_attribute<1>(it, m_radial_map.template create_attribute<1>()); }
    }

    // Assign values
    for (int i=0; i<centroids.size(); ++i)
    {
      auto u=centroids[i];
      bool first_run=true;
      for (Dart_handle it=u; first_run || it!=u; it=m_radial_map.next(m_radial_map.opposite2(it)))
      {
        first_run=false;
        m_radial_map.template info<1>(it)=i;
      }
    }

    // Initialize 0-attributes of m_radial_map
    for (auto it=m_radial_map.darts().begin(), itend=m_radial_map.darts().end(); it!=itend; ++it)
    {
      if (m_radial_map.template attribute<0>(it)==NULL)
      { m_radial_map.template set_attribute<0>(it, m_radial_map.template create_attribute<0>(-1)); }
    }
    for (auto att_it=m_gmap.template attributes<0>().begin(),
         att_itend = m_gmap.template attributes<0>().end(); att_it != att_itend; ++att_it)
    {
      auto it_radial = origin_to_radial[m_copy_to_origin[att_it->dart()]];
      m_radial_map.template info<0>(it_radial)=att_it->info();
    }

    // Remove the marked edges of m_radial_map
    for (auto dh : edge_list)
    {
      CGAL_assertion(m_radial_map.template is_removable<1>(dh));
      m_radial_map.template remove_cell<1>(dh);
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] Facewidth constructor: "<<t.time()<<" seconds."<<std::endl;
    }
  }

  std::vector<Original_dart_const_handle> compute_facewidth(bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_cycle.clear();
    // Find edgewidth of the radial map
    SNC snc_to_find_facewidth(m_radial_map, display_time);
    Path_on_surface<Local_map> edgewidth_of_radial_map=
        snc_to_find_facewidth.compute_edgewidth(display_time);

    int last_vertex_index=-1;
    int last_face_index=-1;
    for (int i=0, n=edgewidth_of_radial_map.length(); i<=n; i++)
    {
      Dart_const_handle dh=edgewidth_of_radial_map[i%n];
      int face_index=m_radial_map.template info<1>(dh);
      if (m_radial_map.template info<0>(dh)==-1) { dh=m_radial_map.next(dh); }
      CGAL_assertion(m_radial_map.template info<0>(dh)!=-1);
      int vertex_index=m_radial_map.template info<0>(dh);

      if (last_face_index==face_index)
      {
        CGAL_assertion(m_radial_map.template belong_to_same_cell<0>(m_vertex_list[last_vertex_index],
                                                                    dh));
        
        m_cycle.push_back(m_copy_to_origin[m_vertex_list[last_vertex_index]]);
        m_cycle.push_back(m_copy_to_origin[m_face_list[face_index]]);
      }

      // if (first_vertex_index == -1) first_vertex_index = vertex_index;
      last_vertex_index=vertex_index;
      last_face_index=face_index;
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_facewidth: "<<t.time()<<" seconds."<<std::endl;
    }

    return m_cycle;
  }

protected:
  Local_map m_radial_map;
  std::vector<Dart_handle> m_vertex_list, m_face_list;
  std::vector<Original_dart_const_handle> m_cycle;
  typename Original_map_wrapper::Origin_to_copy_map m_origin_to_copy;
  typename Original_map_wrapper::Copy_to_origin_map m_copy_to_origin;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
