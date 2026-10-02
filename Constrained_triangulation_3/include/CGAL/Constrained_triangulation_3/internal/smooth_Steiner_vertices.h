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

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
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
  template <typename Cell_handle>
  struct Cells_pmap
  {
    std::shared_ptr<std::unordered_set<Cell_handle>>
      cells_{ std::make_shared<std::unordered_set<Cell_handle>>() };

    using category = boost::read_write_property_map_tag;
    using reference = bool;
    using value_type = bool;
    using key_type = Cell_handle;

    friend reference get(const Cells_pmap& m, key_type c) { return m.cells_->find(c) != m.cells_->end(); }
    friend void put(Cells_pmap& m, key_type c, value_type value)
    {
      if(value)
        m.cells_->insert(c);
      else if(m.cells_->find(c) != m.cells_->end())
        m.cells_->erase(c);
    }
  };

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
      return type == CDT_3_vertex_type::INPUT_VERTEX
          || type == CDT_3_vertex_type::BBOX;
    }
  };

  template<typename CDT>
  struct Is_finite_cell_pmap
  {
    using category = boost::read_write_property_map_tag;
    using reference = bool;
    using value_type = bool;
    using key_type = typename CDT::Cell_handle;

    const typename CDT::Vertex_handle infinite_vertex_;

    friend reference get(const Is_finite_cell_pmap& m, key_type cell)
    {
      for (int i = 0; i < 4; ++i)
      {
        if(cell->vertex(i) == m.infinite_vertex_)
          return false;
      }
      return true;
    }
    friend void put(Is_finite_cell_pmap&, const key_type, const value_type){}
  };

  template <typename Vertex_handle>
  struct Is_not_volume_Steiner_point_map
  {
    using category = boost::readable_property_map_tag;
    using reference = bool;
    using value_type = bool;
    using key_type = Vertex_handle;

    friend reference get(const Is_not_volume_Steiner_point_map&, key_type v)
    {
      auto type = v->ccdt_3_data().vertex_type();
      return type != CDT_3_vertex_type::STEINER_IN_VOLUME;
    }
  };
} // end namespace internal

template <typename CDT_3>
std::size_t refine_around_Steiner_points_in_volume(CDT_3& cdt)
{
  std::vector<typename CDT_3::Vertex_handle> steiner_in_volume;

  for(auto v : cdt.triangulation().finite_vertex_handles())
    if(v->ccdt_3_data().vertex_type() == CDT_3_vertex_type::STEINER_IN_VOLUME)
      steiner_in_volume.push_back(v);

  std::size_t nb_insertions = 0;
  for(auto v : steiner_in_volume)
  {
    std::vector<typename CDT_3::Cell_handle> cells;
    cdt.triangulation().incident_cells(v, std::back_inserter(cells));
    bool do_refine = false;
    for(auto c : cells)
    {
      if(cdt.triangulation().is_infinite(c))
        continue;
      if(CGAL::volume(c->vertex(0)->point(),
                      c->vertex(1)->point(),
                      c->vertex(2)->point(),
                      c->vertex(3)->point()) < 0.)
      {
        do_refine = true;
        break;
      }
    }

    if(do_refine)
    {
      for(auto c : cells)
      {
        if(cdt.triangulation().is_infinite(c))
          continue;
        const auto p = CGAL::centroid(c->vertex(0)->point(),
                                      c->vertex(1)->point(),
                                      c->vertex(2)->point(),
                                      c->vertex(3)->point());
        auto new_v = cdt.triangulation().insert(p);
        new_v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_IN_VOLUME);
        ++nb_insertions;
      }
    }
  }
  return nb_insertions;
}

template <typename CDT_3>
std::size_t refine_cells(CDT_3& cdt)
{
  std::size_t nb_insertions = 0;

  struct Cell_and_point
  {
    typename CDT_3::Cell_handle c;
    typename CDT_3::Triangulation::Geom_traits::Point_3 p;
  };

  std::vector<Cell_and_point> cells;
  for(auto c : cdt.triangulation().finite_cell_handles())
  {
    const auto p = CGAL::centroid(c->vertex(0)->point(),
                                  c->vertex(1)->point(),
                                  c->vertex(2)->point(),
                                  c->vertex(3)->point());
    cells.push_back({c, p});
  }

  struct Facet_from_outside
  {
    bool in_complex;
    typename CDT_3::Facet facet;
    typename CDT_3::Surface_patch_index patch = {};
  };

  for(const auto& [c, p] : cells)
  {
    //save facets seen from outside c
    std::array<Facet_from_outside, 4> facets_from_outside;
    for (int i = 0; i < 4; ++i)
    {
      auto fi = cdt.triangulation().mirror_facet({c, i});
      if(cdt.is_in_complex(fi))
        facets_from_outside[i] = {true, fi, cdt.surface_patch_index(fi)};
      else
        facets_from_outside[i] = {false, fi, typename CDT_3::Surface_patch_index{}};
    }

    //clear indices from c
    for(int i = 0; i < 4; ++i)
      if(cdt.is_in_complex(c, i))
        cdt.remove_from_complex(c, i);

    //insert steiner point
    auto new_v = cdt.triangulation().insert(p);
    new_v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_IN_VOLUME);
    ++nb_insertions;

    //restore patches for facets seen from outside c
    std::vector<typename CDT_3::Cell_handle> new_cells;
    cdt.triangulation().finite_incident_cells(new_v, std::back_inserter(new_cells));
    for(auto nc : new_cells)
    {
      for(int i = 0; i < 4; ++i)
      {
        auto nfi = cdt.triangulation().mirror_facet({nc, i});
        // restore index
        if(facets_from_outside[nfi.second].in_complex)
          cdt.add_to_complex(nfi, facets_from_outside[nfi.second].patch);
      }
    }
  }
  return nb_insertions;
}

template<typename CDT_3>
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

template<typename CDT_3>
void remove_from_complex_and_reset_index(const typename CDT_3::Facet& f, CDT_3& cdt)
{
  if(cdt.is_in_complex(f))
    cdt.remove_from_complex(f);
  cdt.set_surface_patch_index(f, typename CDT_3::Surface_patch_index{});
}

template<typename CDT_3>
auto split_internal_edge(const typename CDT_3::Edge& e, CDT_3& cdt)
{
  using Vertex_handle = typename CDT_3::Vertex_handle;
  using Cell_handle = typename CDT_3::Cell_handle;
  using Facet = typename CDT_3::Facet;
  using Edge = typename CDT_3::Edge;
  using Patch_index = typename CDT_3::Surface_patch_index;

  const Vertex_handle v1 = e.first->vertex(e.second);
  const Vertex_handle v2 = e.first->vertex(e.third);
  Vertex_handle new_v = cdt.triangulation().tds().insert_in_edge(e);
  new_v->set_point(CGAL::midpoint(v1->point(), v2->point()));
  new_v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_IN_VOLUME);

  // restore facets seen from outside the cells incident to the edge
  std::vector<Cell_handle> new_cells;
  cdt.triangulation().finite_incident_cells(new_v, std::back_inserter(new_cells));
  for(Cell_handle new_c : new_cells)
  {
    const int index_new_v = new_c->index(new_v);
    for(int i = 0; i < 4; ++i)
    {
      const Facet fi(new_c, i);
      const Facet mfi = cdt.triangulation().mirror_facet(fi);

      const Patch_index patch = cdt.surface_patch_index(mfi); // backup patch index

      remove_from_complex_and_reset_index(fi, cdt);
      remove_from_complex_and_reset_index(mfi, cdt);

      // outer hull facets
      if(i == index_new_v && patch != Patch_index{})
        cdt.add_to_complex(mfi, patch);
    }
  }
  return new_v;
}

template <typename CDT_3>
std::size_t refine_internal_edges(CDT_3& cdt)
{
  using Edge = typename CDT_3::Edge;

  std::size_t nb_insertions = 0;
  struct Edge_vv
  {
    typename CDT_3::Vertex_handle v1, v2;
  };
  auto make_edge_vv
    = [](typename CDT_3::Vertex_handle v1,
         typename CDT_3::Vertex_handle v2)
          {
            if(v1 < v2)
              return Edge_vv{v1, v2};
            else
              return Edge_vv{v2, v1};
          };

  std::vector<Edge_vv> edges;
  for(auto e : cdt.triangulation().finite_edges())
  {
    if(cdt.is_in_complex(e))
      continue;

    auto v1 = e.first->vertex(e.second);
    auto v2 = e.first->vertex(e.third);
    if(v1->ccdt_3_data().vertex_type() == CDT_3_vertex_type::BBOX ||
       v2->ccdt_3_data().vertex_type() == CDT_3_vertex_type::BBOX)
      continue;

    edges.push_back(make_edge_vv(v1, v2));
  }

  std::ofstream ofs("edges_to_refine.polylines.txt");
  for(const auto& [v1, v2] : edges)
    ofs << "2 " << v1->point() << " " << v2->point() << std::endl;
  ofs.close();

  std::ofstream ofs2("edges_to_refine_steiner_points.xyz");
  ofs2.precision(17);
  CGAL::Static_boolean_property_map<typename CDT_3::Cell_handle, true> select;
  for(const auto& [v1, v2] : edges)
  {
    typename CDT_3::Cell_handle c;
    int i, j;
    if (cdt.triangulation().tds().is_edge(v1, v2, c, i, j))
    {
      //auto new_v = CGAL::Tetrahedral_remeshing::internal::split_edge({c, i, j}, select, cdt);
      auto new_v = split_internal_edge(Edge(c, i, j), cdt);
      //if(new_v != typename CDT_3::Vertex_handle())
      {
        //new_v->ccdt_3_data().set_vertex_type(CDT_3_vertex_type::STEINER_IN_VOLUME);
        ++nb_insertions;
        ofs2 << new_v->point() << std::endl;
      }
    }

    std::ofstream ofs00("out_after_edge_split.mesh");
    ofs00.precision(17);
    CGAL::IO::write_MEDIT(ofs00, cdt);
    ofs00.close();
    std::cout << "Refined edge " << v1->point() << " " << v2->point() << std::endl;
    std::cout << "File " << "out_after_edge_split.mesh" << " written." << std::endl;
    std::cout << std::endl;
  }
  ofs2.close();

  return nb_insertions;
}

template <typename CDT_3>
std::size_t perturb_slivers(CDT_3& cdt, const double bbox_max_span)
{
  const double min_volume = bbox_max_span * bbox_max_span * bbox_max_span * 1e-8;
  std::cout << "Perturb slivers with volume < " << min_volume << std::endl;

  using GT = typename CDT_3::Triangulation::Geom_traits;
  using Vector_3 = typename GT::Vector_3;

  CGAL::Random rng;
  std::size_t nb_perturbations = 0;
  for(auto v : cdt.triangulation().finite_vertex_handles())
  {
    if(v->ccdt_3_data().vertex_type() != CDT_3_vertex_type::STEINER_IN_VOLUME)
      continue;
    std::vector<typename CDT_3::Cell_handle> cells;
    cdt.triangulation().finite_incident_cells(v, std::back_inserter(cells));
    bool do_perturb = false;
    for(auto c : cells)
    {
      if(CGAL::abs(CGAL::volume(c->vertex(0)->point(),
                                c->vertex(1)->point(),
                                c->vertex(2)->point(),
                                c->vertex(3)->point())) < min_volume)
      {
        do_perturb = true;
        break;
      }
    }
    if(do_perturb)
    {
      const auto& p = v->point();
      const Vector_3 perturbation(rng.get_double(-1e-2, 1e-2) * p.x(),
                                  rng.get_double(-1e-2, 1e-2) * p.y(),
                                  rng.get_double(-1e-2, 1e-2) * p.z());
      const auto p_perturbed = p + perturbation;
      v->set_point(p_perturbed);
      ++nb_perturbations;
    }
  }
  return nb_perturbations;
}

template <typename CDT_3>
void dump_negative_tetrahedra(const char* filename, const CDT_3& cdt)
{
  std::ofstream ofs(filename);
  ofs.precision(17);
  for(auto c : cdt.finite_cell_handles()) {
    if(CGAL::volume(c->vertex(0)->point(),
                    c->vertex(1)->point(),
                    c->vertex(2)->point(),
                    c->vertex(3)->point()) < 0.)
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
                                         const bool remesh,
                                         const bool refine_and_smooth)
  {
    using Cell_handle = typename CDT_3::Cell_handle;
    using Vertex_handle = typename CDT_3::Vertex_handle;

    std::cout << "#Negative tetrahedra before Steiner smoothing: "
              << count_negative_tetrahedra(cdt) << std::endl;

    std::ofstream ofs00("out_before_steiner_smoothing.mesh");
    ofs00.precision(17);
    CGAL::IO::write_MEDIT(ofs00, cdt, CGAL::parameters::all_cells(true));
    ofs00.close();

    if (count_negative_tetrahedra(cdt) == 0)
    {
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
      if(count_negative_tetrahedra(cdt) == 0)
      {
        std::cout << "Good :-) No negative tetrahedra left after collapse!" << std::endl;
        return;
      }

      // split long internal edges
      const double target_edge_length = 0.0005 * bbox_max_span;
      CGAL::tetrahedral_isotropic_remeshing(cdt, target_edge_length,
                                            CGAL::parameters::number_of_iterations(5)
                                                .nb_smoothing_iterations(5)
                                                .nb_flip_smooth_iterations(5)
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

    //if(refine_and_smooth)
    for(int i = 0; i < 2; ++i)
    {
      std::cout << "Refining around Steiner points in volume...";
      std::cout.flush();
      std::size_t nb_insertions = refine_cells(cdt);

      std::ofstream ofs0("out_after_cells_refinement.mesh");
      ofs0.precision(17);
      CGAL::IO::write_MEDIT(ofs0, cdt);
      ofs0.close();
      dump_negative_tetrahedra("out_negative_tetrahedra_after_cells_refinement.polylines.txt", cdt);

//      nb_insertions += refine_internal_edges(cdt);
        //refine_around_Steiner_points_in_volume(cdt);
      std::cout << "done (" << nb_insertions << " points inserted)" << std::endl;
      std::cout << cdt.statistics() << std::endl;

      std::ofstream ofs("out_after_refinement.mesh");
      ofs.precision(17);
      CGAL::IO::write_MEDIT(ofs, cdt);
      ofs.close();
      dump_negative_tetrahedra("out_negative_tetrahedra_after_refinement.polylines.txt", cdt);

//      std::cout << std::endl;
//      std::cout << "Perturb slivers...";
//      std::cout.flush();
//      std::size_t nb_perturbations = perturb_slivers(cdt, bbox_max_span);
//      std::cout << "done (" << nb_perturbations << " perturbations)" << std::endl;
//      std::cout << cdt.statistics() << std::endl;
//
//      std::ofstream ofs1("out_after_perturb.mesh");
//      ofs1.precision(17);
//      CGAL::IO::write_MEDIT(ofs1, cdt);
//      ofs1.close();
//      dump_negative_tetrahedra
//      ("out_negative_tetrahedra_after_perturb.polylines.txt", cdt);

      std::cout << std::endl;
      std::cout << "Mesh smoothing around Steiner points in volume...";
      std::cout.flush();
      internal::Is_not_volume_Steiner_point_map<Vertex_handle> is_not_volume_Steiner_point_pmap;
      CGAL::boundary_aware_mesh_smoothing(
          cdt,
          CGAL::Mesh_smoothing_3::C3t3_no_projection<CDT_3>(),
          CGAL::parameters::vertex_is_constrained_map(is_not_volume_Steiner_point_pmap)
          .number_of_iterations(1500));
      std::cout << "done." << std::endl;
      std::cout << cdt.statistics() << std::endl;

      std::ofstream ofs2("out_after_refinement_and_smoothing.mesh");
      ofs2.precision(17);
      CGAL::IO::write_MEDIT(ofs2, cdt);
      ofs2.close();
      dump_negative_tetrahedra("out_negative_tetrahedra_after_refinement_and_smoothing.polylines.txt", cdt);

      if(count_negative_tetrahedra(cdt) == 0)
      {
        std::cout << "Good :-) No negative tetrahedra left after smoothing!" << std::endl;
        return;
      }
    }

    std::cout << cdt.statistics() << std::endl;
  }
} // namespace CGAL

#endif // CGAL_INTERNAL_CDT_3_SMOOTH_STEINER_VERTICES_H
