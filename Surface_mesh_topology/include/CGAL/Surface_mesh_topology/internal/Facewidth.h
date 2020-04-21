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

  using Original_map_wrapper=internal::Generic_map_selector<Mesh, Items_for_shortest_noncontractible_cycle>;
  using Original_dart_const_handle=typename Original_map_wrapper::Dart_const_handle_original;

  using Local_map        =typename Original_map_wrapper::Generic_map;
  using Dart_handle      =typename Local_map::Dart_handle;
  using Dart_const_handle=typename Local_map::Dart_const_handle;
  using size_type        = typename Local_map::size_type;

  using Path             = CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;
  using SNC              = Shortest_noncontractible_cycle<Local_map, false>;

  Facewidth(const Mesh& amesh, bool display_time=false):
    m_snc_to_find_facewidth(nullptr)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    typename Original_map_wrapper::Origin_to_copy_map origin_to_radial;
    typename Original_map_wrapper::Copy_to_origin_map local_radial_to_origin;

    m_is_perforated=m_radial_map.get_new_mark();
    Original_map_wrapper::copy(m_radial_map, amesh, origin_to_radial,
                               local_radial_to_origin, m_is_perforated);

    m_radial_map.negate_mark(m_is_perforated);
    // Remove all boundary by adding faces, marked with m_is_perforated
    m_radial_map.template close<2>();
    m_radial_map.negate_mark(m_is_perforated);

    m_mark_new_vertex=m_radial_map.get_new_mark();
    size_type face_treated=m_radial_map.get_new_mark();
    size_type original_darts=m_radial_map.get_new_mark();
    m_radial_map.negate_mark(original_darts); // All original darts are marked

    for (typename Local_map::Dart_range::iterator it=m_radial_map.darts().begin();
         it!=m_radial_map.darts().end(); ++it)
    {
      if (!m_radial_map.is_marked(it, m_is_perforated) &&
          !m_radial_map.is_marked(it, face_treated) &&
          m_radial_map.is_marked(it, original_darts))
      {
        m_radial_map.template mark_cell<2>(it, face_treated);
        Dart_handle new_vertex=m_radial_map.insert_cell_0_in_cell_2(it);
        for (auto itv=m_radial_map.template darts_of_cell_basic<0>(new_vertex, m_mark_new_vertex).begin(),
             itvend=m_radial_map.template darts_of_cell_basic<0>(new_vertex, m_mark_new_vertex).end();
             itv!=itvend; ++itv)
        {
          m_radial_map.mark(itv, m_mark_new_vertex);
          m_radial_to_original[m_radial_map.opposite2(itv)]=
              local_radial_to_origin[m_radial_map.next(itv)];
        }
      }
    }

    // Remove the original edges of m_radial_map, and create vertex info.
    for (auto it=m_radial_map.darts().begin(), itend=m_radial_map.darts().end();
         it!=itend; ++it)
    {
      if (m_radial_map.template attribute<0>(it)==nullptr)
      { m_radial_map.template set_attribute<0>
          (it, m_radial_map.template create_attribute<0>(-1)); }

      if (m_radial_map.is_marked(it, original_darts))
      { m_radial_map.template remove_cell<1>(it); }
    }

    CGAL_assertion(m_radial_map.is_whole_map_unmarked(original_darts));
    CGAL_assertion(m_radial_map.is_whole_map_unmarked(face_treated));

    m_radial_map.free_mark(original_darts);
    m_radial_map.free_mark(face_treated);

    m_snc_to_find_facewidth=std::make_unique<SNC>(&m_radial_map, m_is_perforated);

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] Facewidth constructor: "<<t.time()<<" seconds."<<std::endl;
    }
  }

  std::vector<Original_dart_const_handle> compute_face_width(bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    // Find edgewidth of the radial map
    Path_on_surface<Local_map> edgewidth_of_radial_map=
        m_snc_to_find_facewidth->compute_edge_width();

    std::vector<Original_dart_const_handle> cycle;
    cycle.reserve(edgewidth_of_radial_map.length());
    for (std::size_t i=0, n=edgewidth_of_radial_map.length(); i<n; i++)
    {
      Dart_handle dh=m_radial_map.dart_handle
          (const_cast<typename Local_map::Dart&>(*edgewidth_of_radial_map.get_ith_real_dart(i)));
      if (!m_radial_map.is_marked(dh, m_mark_new_vertex))
      { cycle.push_back(m_radial_to_original[dh]); }
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_face_width: "<<t.time()<<" seconds."<<std::endl;
    }

    return cycle;
  }

protected:
  Local_map m_radial_map;
  size_type m_is_perforated;   /// mark for perforated darts
  size_type m_mark_new_vertex; /// mark for new vertices
  typename Original_map_wrapper::Copy_to_origin_map m_radial_to_original;
  std::unique_ptr<SNC> m_snc_to_find_facewidth;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
