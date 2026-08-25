// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include "test_common.h"

#include <algorithm>
#include <cassert>

using namespace ms3_test;

namespace {

using Projection_mode = CGAL::Mesh_smoothing_3::Projection_weight_mode;

template <typename C3t3T>
bool has_query(Recording_cts<C3t3T> const& cts,
               typename Recording_cts<C3t3T>::Patch_face const& q)
{
  return std::find(cts.patch_queries.begin(), cts.patch_queries.end(), q) != cts.patch_queries.end();
}

template <typename C3t3T>
bool has_query(Recording_cts<C3t3T> const& cts,
               typename Recording_cts<C3t3T>::Curve_edge const& q)
{
  return std::find(cts.curve_queries.begin(), cts.curve_queries.end(), q) != cts.curve_queries.end();
}

template <typename C3t3T>
void expect_only_expected_queries(
    Recording_cts<C3t3T> const& cts,
    std::vector<typename Recording_cts<C3t3T>::Patch_face> const& expected_patches,
    std::vector<typename Recording_cts<C3t3T>::Curve_edge> const& expected_curves)
{
  for (auto const& q : cts.patch_queries) {
    assert(std::find(expected_patches.begin(), expected_patches.end(), q) != expected_patches.end());
  }
  for (auto const& q : cts.curve_queries) {
    assert(std::find(expected_curves.begin(), expected_curves.end(), q) != expected_curves.end());
  }
  for (auto const& expected : expected_patches) {
    assert(has_query(cts, expected));
  }
  for (auto const& expected : expected_curves) {
    assert(has_query(cts, expected));
  }
}

template <typename C3t3T>
std::vector<typename C3t3T::Cell_handle> collect_cells_in_complex(C3t3T const& c3t3)
{
  return std::vector<typename C3t3T::Cell_handle>(c3t3.cells_in_complex().begin(),
                                                  c3t3.cells_in_complex().end());
}

template <typename C3t3T>
std::vector<typename C3t3T::Facet> collect_facets_in_complex(C3t3T const& c3t3)
{
  return std::vector<typename C3t3T::Facet>(c3t3.facets_in_complex().begin(),
                                            c3t3.facets_in_complex().end());
}

template <typename C3t3T>
std::vector<typename C3t3T::Edge> collect_edges_in_complex(C3t3T const& c3t3)
{
  return std::vector<typename C3t3T::Edge>(c3t3.edges_in_complex().begin(),
                                           c3t3.edges_in_complex().end());
}

template <typename C3t3T>
std::vector<typename C3t3T::Vertex_handle> collect_vertices_in_complex(C3t3T const& c3t3)
{
  return std::vector<typename C3t3T::Vertex_handle>(c3t3.vertices_in_complex().begin(),
                                                    c3t3.vertices_in_complex().end());
}

template <typename C3t3T>
void assert_structure_preserved(C3t3T const& before, C3t3T const& after)
{
  assert(before.triangulation().number_of_vertices() == after.triangulation().number_of_vertices());
  assert(before.triangulation().number_of_finite_cells() == after.triangulation().number_of_finite_cells());
  assert(before.c3t3.number_of_cells_in_complex() == after.c3t3.number_of_cells_in_complex());
  assert(before.c3t3.number_of_facets_in_complex() == after.c3t3.number_of_facets_in_complex());
  assert(before.c3t3.number_of_edges_in_complex() == after.c3t3.number_of_edges_in_complex());
  assert(before.c3t3.number_of_vertices_in_complex() == after.c3t3.number_of_vertices_in_complex());

  auto before_cells = collect_cells_in_complex(before.c3t3);
  auto after_cells = collect_cells_in_complex(after.c3t3);
  assert(before_cells.size() == after_cells.size());
  for (auto c : before_cells) {
    assert(has_value(after_cells, c));
    assert(c->subdomain_index() == after.c3t3.subdomain_index(c));
    for (int i = 0; i < 4; ++i) {
      assert(c->vertex(i) == c->vertex(i));
    }
  }

  auto before_facets = collect_facets_in_complex(before.c3t3);
  auto after_facets = collect_facets_in_complex(after.c3t3);
  assert(before_facets.size() == after_facets.size());
  for (auto const& f : before_facets) {
    assert(has_value(after_facets, f));
    assert(before.c3t3.surface_patch_index(f) == after.c3t3.surface_patch_index(f));
  }

  auto before_edges = collect_edges_in_complex(before.c3t3);
  auto after_edges = collect_edges_in_complex(after.c3t3);
  assert(before_edges.size() == after_edges.size());
  for (auto const& e : before_edges) {
    assert(has_value(after_edges, e));
    assert(before.c3t3.curve_index(e) == after.c3t3.curve_index(e));
  }

  auto before_vertices = collect_vertices_in_complex(before.c3t3);
  auto after_vertices = collect_vertices_in_complex(after.c3t3);
  assert(before_vertices.size() == after_vertices.size());
  for (auto v : before_vertices) {
    assert(has_value(after_vertices, v));
    assert(before.c3t3.corner_index(v) == after.c3t3.corner_index(v));
    assert(before.c3t3.in_dimension(v) == after.c3t3.in_dimension(v));
  }
}

template <typename C3t3T>
void install_default_queries(Recording_cts<C3t3T>& cts)
{
  using GT = typename Recording_cts<C3t3T>::Geom_traits;
  using Vector_3 = typename GT::Vector_3;
  using Tangent_space = typename Recording_cts<C3t3T>::Tangent_space;

  cts.patch_results.emplace(11, Tangent_space{Point_3(0., 0., 0.), Vector_3(0., 0., 1.),
                                              Projection_mode::DEFAULT, 1.});
  cts.patch_results.emplace(12, Tangent_space{Point_3(0., 0., 0.), Vector_3(0., 0., 1.),
                                              Projection_mode::STRONG, 1.});
  cts.curve_results.emplace(21, Tangent_space{Point_3(1.3, 0.1, 0.), Vector_3(1., 0., 0.),
                                              Projection_mode::SOFT, 1.});
  cts.curve_results.emplace(22, Tangent_space{Point_3(0.2, 1.1, 0.), Vector_3(0., 1., 0.),
                                              Projection_mode::CUSTOM, 0.0});
}

void test_classification_and_structure()
{
  auto fixture = make_fixture(false, false);
  auto before_cells = collect_cells_in_complex(fixture.c3t3);
  auto before_facets = collect_facets_in_complex(fixture.c3t3);
  auto before_edges = collect_edges_in_complex(fixture.c3t3);
  auto before_vertices = collect_vertices_in_complex(fixture.c3t3);
  std::map<Cell_handle, Subdomain_index> before_subdomains;
  std::map<Facet, Surface_patch_index> before_patches;
  std::map<Edge, Curve_index> before_curves;
  std::map<Vertex_handle, Corner_index> before_corners;
  for (auto c : before_cells) before_subdomains[c] = fixture.c3t3.subdomain_index(c);
  for (auto const& f : before_facets) before_patches[f] = fixture.c3t3.surface_patch_index(f);
  for (auto const& e : before_edges) before_curves[e] = fixture.c3t3.curve_index(e);
  for (auto v : before_vertices) before_corners[v] = fixture.c3t3.corner_index(v);

  Recording_cts<C3t3> cts;
  install_default_queries(cts);

  auto status = CGAL::boundary_aware_mesh_smoothing(
      fixture.c3t3, cts, CGAL::parameters::number_of_iterations(1));

  assert_all_finite(fixture.c3t3);
  assert(status.valid_mesh() == (status.nb_invalid_elements == 0));

  assert(before_cells.size() == collect_cells_in_complex(fixture.c3t3).size());
  assert(before_facets.size() == collect_facets_in_complex(fixture.c3t3).size());
  assert(before_edges.size() == collect_edges_in_complex(fixture.c3t3).size());
  assert(before_vertices.size() == collect_vertices_in_complex(fixture.c3t3).size());
  for (auto c : before_cells) {
    assert(before_subdomains[c] == fixture.c3t3.subdomain_index(c));
  }
  for (auto const& f : before_facets) {
    assert(before_patches[f] == fixture.c3t3.surface_patch_index(f));
  }
  for (auto const& e : before_edges) {
    assert(before_curves[e] == fixture.c3t3.curve_index(e));
  }
  for (auto v : before_vertices) {
    assert(before_corners[v] == fixture.c3t3.corner_index(v));
  }

  expect_only_expected_queries(
      cts,
      {std::make_pair(11u, fixture.shared_facet), std::make_pair(12u, fixture.complex_facet_2)},
      {std::make_pair(21u, fixture.curve_edge_1), std::make_pair(22u, fixture.curve_edge_2)});

  assert(!has_query(cts, std::make_pair(11u, fixture.unmarked_other_facet)));

  auto edge_ab = fixture.unmarked_edge;
  auto edge_ac = find_edge(fixture.c3t3.triangulation(), fixture.v[0], fixture.v[2]);
  auto edge_bc = find_edge(fixture.c3t3.triangulation(), fixture.v[1], fixture.v[2]);
  assert(!has_query(cts, std::make_pair(21u, edge_ab)));
  assert(!has_query(cts, std::make_pair(21u, edge_ac)));
  assert(!has_query(cts, std::make_pair(21u, edge_bc)));
}

void test_zero_bounds_and_invalid_counting()
{
  {
    auto fixture = make_fixture(false, false);
    auto before = finite_vertices(fixture.c3t3);
    Recording_cts<C3t3> cts;
    install_default_queries(cts);

    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::number_of_iterations(1).max_number_of_evaluations(0));

    assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::MAX_NUMBER_OF_METRIC_EVALUATIONS_REACHED);
    assert(status.nb_iterations > 0);
    assert(status.nb_vertex_updates > 0);
    auto after = finite_vertices(fixture.c3t3);
    assert(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
      assert(finite_point(after[i]));
    }
  }

  {
    auto fixture = make_fixture(false, false);
    auto before = finite_vertices(fixture.c3t3);
    Recording_cts<C3t3> cts;
    install_default_queries(cts);

    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::number_of_iterations(1).maximum_running_time(0));

    assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::TIME_LIMIT_REACHED);
    assert(status.nb_iterations > 0);
    assert(status.nb_vertex_updates > 0);
    auto after = finite_vertices(fixture.c3t3);
    assert(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
      assert(finite_point(after[i]));
    }
  }

  {
    auto fixture = make_fixture(false, true);
    auto o = CGAL::orientation(fixture.other_cell->vertex(0)->point(),
                               fixture.other_cell->vertex(1)->point(),
                               fixture.other_cell->vertex(2)->point(),
                               fixture.other_cell->vertex(3)->point());
    assert(o == CGAL::NEGATIVE);

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts, CGAL::parameters::number_of_iterations(0));

    assert(status.nb_initial_invalid_elements == 0);
    assert(status.nb_invalid_elements == 0);
    assert(status.valid_mesh());
    assert(finite_vertices(fixture.c3t3).size() == 5);
  }
}

void test_all_vertices_frozen_status()
{
  auto fixture = make_fixture(false, false);
  std::map<Vertex_handle, bool> vmap;
  for (auto v : fixture.c3t3.triangulation().finite_vertex_handles()) {
    vmap[v] = true;
  }

  Recording_cts<C3t3> cts;
  install_default_queries(cts);
  auto before = finite_vertices(fixture.c3t3);
  auto status = CGAL::boundary_aware_mesh_smoothing(
      fixture.c3t3, cts,
      CGAL::parameters::vertex_is_constrained_map(boost::make_assoc_property_map(vmap))
          .number_of_iterations(5));

  assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::ALL_VERTICES_FROZEN);
  assert(status.nb_vertex_updates == 0);
  auto after = finite_vertices(fixture.c3t3);
  assert(before.size() == after.size());
  for (std::size_t i = 0; i < before.size(); ++i) {
    assert(same_point(before[i], after[i]));
  }
}

void test_hard_constraints_and_frozen_status()
{
  {
    auto fixture = make_fixture(false, false);
    std::map<Vertex_handle, bool> vmap;
    vmap[fixture.v[1]] = true;

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto before_b = fixture.v[1]->point();
    auto before_d = fixture.v[3]->point();
    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::vertex_is_constrained_map(boost::make_assoc_property_map(vmap))
            .number_of_iterations(3));

    assert(same_point(before_b, fixture.v[1]->point()));
    assert(status.nb_vertex_updates > 0);
    assert(!same_point(before_d, fixture.v[3]->point()));
  }

  {
    auto fixture = make_fixture(false, false);
    std::map<std::pair<Vertex_handle, Vertex_handle>, bool> emap;
    emap[std::make_pair(fixture.v[3], fixture.v[1])] = true;

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto before_b = fixture.v[1]->point();
    auto before_d = fixture.v[3]->point();
    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::edge_is_constrained_map(boost::make_assoc_property_map(emap))
            .number_of_iterations(1));

    (void)status;
    assert(same_point(before_b, fixture.v[1]->point()));
    assert(same_point(before_d, fixture.v[3]->point()));
  }

  {
    auto fixture = make_fixture(false, false);
    std::map<std::pair<Vertex_handle, Vertex_handle>, bool> emap;
    emap[std::make_pair(fixture.v[1], fixture.v[3])] = true;

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto before_b = fixture.v[1]->point();
    auto before_d = fixture.v[3]->point();
    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::edge_is_constrained_map(boost::make_assoc_property_map(emap))
            .number_of_iterations(1));

    (void)status;
    assert(same_point(before_b, fixture.v[1]->point()));
    assert(same_point(before_d, fixture.v[3]->point()));
  }

  {
    auto fixture = make_fixture(false, false);
    std::map<Facet, bool> fmap;
    fmap[fixture.c3t3.triangulation().mirror_facet(fixture.shared_facet)] = true;

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto before_a = fixture.v[0]->point();
    auto before_b = fixture.v[1]->point();
    auto before_c = fixture.v[2]->point();

    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::facet_is_constrained_map(boost::make_assoc_property_map(fmap))
            .number_of_iterations(1));

    assert(status.nb_vertex_updates > 0);
    assert(same_point(before_a, fixture.v[0]->point()));
    assert(same_point(before_b, fixture.v[1]->point()));
    assert(same_point(before_c, fixture.v[2]->point()));
  }

  {
    auto fixture = make_fixture(false, false);
    std::map<Vertex_handle, bool> vmap;
    std::map<std::pair<Vertex_handle, Vertex_handle>, bool> emap;
    std::map<Facet, bool> fmap;

    vmap[fixture.v[1]] = true;
    emap[std::make_pair(fixture.v[2], fixture.v[3])] = true;
    fmap[fixture.shared_facet] = true;

    Recording_cts<C3t3> cts;
    install_default_queries(cts);
    auto before = finite_vertices(fixture.c3t3);
    auto status = CGAL::boundary_aware_mesh_smoothing(
        fixture.c3t3, cts,
        CGAL::parameters::vertex_is_constrained_map(boost::make_assoc_property_map(vmap))
            .edge_is_constrained_map(boost::make_assoc_property_map(emap))
            .facet_is_constrained_map(boost::make_assoc_property_map(fmap))
            .number_of_iterations(5));

    assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::ALL_VERTICES_FROZEN);
    assert(status.nb_vertex_updates == 0);
    auto after = finite_vertices(fixture.c3t3);
    assert(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
      assert(same_point(before[i], after[i]));
    }
  }
}

} // namespace

int main()
{
  test_classification_and_structure();
  test_zero_bounds_and_invalid_counting();
  test_all_vertices_frozen_status();
  test_hard_constraints_and_frozen_status();
  return EXIT_SUCCESS;
}
