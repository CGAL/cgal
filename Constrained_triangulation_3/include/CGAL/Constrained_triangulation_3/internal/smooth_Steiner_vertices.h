// Copyright (c) 2019-2026  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#ifndef CGAL_INTERNAL_CDT_3_SMOOTH_STEINER_VERTICES_H
#define CGAL_INTERNAL_CDT_3_SMOOTH_STEINER_VERTICES_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>
#include <CGAL/Random.h>

#include <CGAL/IO/write_MEDIT.h>

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <memory>

namespace CGAL
{
namespace internal
{

template <typename Vertex_handle>
struct Is_input_vertex_map
{
  using category = boost::readable_property_map_tag;
  using reference = bool;
  using value_type = bool;
  using key_type = Vertex_handle;
  friend reference get(const Is_input_vertex_map&, key_type v)
  {
    auto type = v->ccdt_3_data().vertex_type();
    return type == CDT_3_vertex_type::INPUT_VERTEX || type == CDT_3_vertex_type::BBOX;
  }
};

template <typename CDT>
struct Is_finite_cell_pmap
{
  using category = boost::read_write_property_map_tag;
  using reference = bool;
  using value_type = bool;
  using key_type = typename CDT::Cell_handle;

  const typename CDT::Vertex_handle infinite_vertex_;

  friend reference get(const Is_finite_cell_pmap& m, key_type cell)
  {
    for(int i = 0; i < 4; ++i)
    {
      if(cell->vertex(i) == m.infinite_vertex_)
        return false;
    }
    return true;
  }
  friend void put(Is_finite_cell_pmap&, const key_type, const value_type) {}
};

template <typename CDT_3>
std::size_t count_negative_tetrahedra(const CDT_3& cdt)
{
  const auto& tr = cdt.triangulation();
  return std::count_if(
           tr.finite_cells_begin(),
           tr.finite_cells_end(),
           [&](const auto& c)
           {
             return CGAL::POSITIVE != CGAL::orientation(tr.point(c.vertex(0)),
                                                        tr.point(c.vertex(1)),
                                                        tr.point(c.vertex(2)),
                                                        tr.point(c.vertex(3)));
           });
}

template <typename CDT_3>
void dump_negative_tetrahedra(const char* filename, const CDT_3& cdt)
{
  std::ofstream ofs(filename);
  ofs.precision(17);
  for(auto c : cdt.finite_cell_handles())
  {
    if(CGAL::volume(c->vertex(0)->point(), c->vertex(1)->point(),
                    c->vertex(2)->point(), c->vertex(3)->point()) < 0.)
    {
      ofs << "2 " << c->vertex(0)->point() << " " << c->vertex(1)->point() << std::endl;
      ofs << "2 " << c->vertex(0)->point() << " " << c->vertex(2)->point() << std::endl;
      ofs << "2 " << c->vertex(0)->point() << " " << c->vertex(3)->point() << std::endl;
      ofs << "2 " << c->vertex(1)->point() << " " << c->vertex(2)->point() << std::endl;
      ofs << "2 " << c->vertex(1)->point() << " " << c->vertex(3)->point() << std::endl;
      ofs << "2 " << c->vertex(2)->point() << " " << c->vertex(3)->point() << std::endl;
    }
  }
  ofs.close();
}

template <typename CDT_3>
void smooth_Steiner_vertices_in_volume(CDT_3& cdt,
                                       const double bbox_max_span,
                                       const bool remesh)
{
  using Cell_handle = typename CDT_3::Cell_handle;
  using Vertex_handle = typename CDT_3::Vertex_handle;

  std::cout << "#Negative tetrahedra before Steiner smoothing: " << count_negative_tetrahedra(cdt) << std::endl;

  std::ofstream ofs00("out_before_steiner_smoothing.mesh");
  ofs00.precision(17);
  CGAL::IO::write_MEDIT(ofs00, cdt, CGAL::parameters::all_cells(true).with_plc_face_id(true));
  ofs00.close();

  if(count_negative_tetrahedra(cdt) == 0) {
    std::cout << "Good :-) No negative tetrahedra left!" << std::endl;
    return;
  }

  cdt.triangulation().may_have_badly_oriented_cells(true);

  if(remesh)
  {
    internal::Is_input_vertex_map<Vertex_handle> is_input_vertex_pmap;
    internal::Is_finite_cell_pmap<CDT_3> is_finite_cell{cdt.triangulation().infinite_vertex()};

    // collapse and smooth only
    CGAL::tetrahedral_isotropic_remeshing(cdt, bbox_max_span,
                                          CGAL::parameters::number_of_iterations(5)
                                              .nb_smoothing_iterations(3)
                                              .nb_flip_smooth_iterations(0)
                                              .do_split(false)
                                              .do_collapse(true)
                                              .do_flip(true)
                                              .remesh_boundaries(false)
                                              .cell_is_selected_map(is_finite_cell)
                                              .vertex_is_constrained_map(is_input_vertex_pmap));

    std::cout << std::endl;
    std::cout << "AFTER COLLAPSE" << std::endl;
    std::cout << cdt.statistics() << std::endl;

    std::ofstream ofs("out_after_collapse.mesh");
    ofs.precision(17);
    CGAL::IO::write_MEDIT(ofs, cdt);
    ofs.close();
    if(count_negative_tetrahedra(cdt) == 0) {
      std::cout << "Good :-) No negative tetrahedra left after collapse!" << std::endl;
      return;
    }

    // split long internal edges
    const double target_edge_length = 0.0005 * bbox_max_span;
    CGAL::tetrahedral_isotropic_remeshing(cdt, target_edge_length,
                                          CGAL::parameters::number_of_iterations(5)
                                              .nb_smoothing_iterations(5)
                                              .nb_flip_smooth_iterations(10)
                                              .do_split(true)
                                              .do_collapse(true)
                                              .do_flip(true)
                                              .remesh_boundaries(false)
                                              .vertex_is_constrained_map(is_input_vertex_pmap));

    std::cout << std::endl;
    std::cout << "AFTER SPLIT AND BEFORE SMOOTH" << std::endl;
    std::cout << cdt.statistics() << std::endl;

    if(count_negative_tetrahedra(cdt) == 0)
    {
      std::cout << "Good :-) No negative tetrahedra left after remeshing!" << std::endl;
      return;
    }
  }
  std::cout << cdt.statistics() << std::endl;
}

} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_CDT_3_SMOOTH_STEINER_VERTICES_H
