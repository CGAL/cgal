// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include "test_common.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/File_medit.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <random>

using namespace ms3_test;

namespace {

constexpr unsigned max_iterations = 100u;
constexpr unsigned untangle_iterations = 100u;
constexpr double untangle_amplitude = 0.05;

C3t3 load_example_mesh()
{
  C3t3 c3t3;
  std::ifstream is(MS3_EXAMPLE_MESH_FILE, std::ios::in);
  assert(is);
  const bool read_ok = CGAL::IO::read_MEDIT(is, c3t3.triangulation());
  assert(read_ok);
  c3t3.rescan_after_load_of_triangulation();
  return c3t3;
}


template <typename C3t3T>
std::vector<Point_3> sorted_vertices(C3t3T const& c3t3)
{
  auto vertices = finite_vertices(c3t3);
  std::sort(vertices.begin(), vertices.end(), [](Point_3 const& a, Point_3 const& b) {
    return std::array<double, 3>{CGAL::to_double(a.x()), CGAL::to_double(a.y()), CGAL::to_double(a.z())} <
           std::array<double, 3>{CGAL::to_double(b.x()), CGAL::to_double(b.y()), CGAL::to_double(b.z())};
  });
  return vertices;
}

template <typename C3t3T>
void assert_same_vertices(C3t3T const& lhs, C3t3T const& rhs)
{
  auto left = sorted_vertices(lhs);
  auto right = sorted_vertices(rhs);
  assert(left.size() == right.size());
  for (std::size_t i = 0; i < left.size(); ++i) {
    assert(same_point(left[i], right[i]));
  }
}

template <typename C3t3T>
void assert_structure_counts_preserved(C3t3T const& before, C3t3T const& after)
{
  assert(before.triangulation().number_of_vertices() == after.triangulation().number_of_vertices());
  assert(before.triangulation().number_of_finite_cells() == after.triangulation().number_of_finite_cells());
  assert(before.number_of_cells_in_complex() == after.number_of_cells_in_complex());
  assert(before.number_of_facets_in_complex() == after.number_of_facets_in_complex());
  assert(before.number_of_edges_in_complex() == after.number_of_edges_in_complex());
  assert(before.number_of_vertices_in_complex() == after.number_of_vertices_in_complex());
}

template <typename C3t3T>
bool any_vertex_moved(C3t3T const& before, C3t3T const& after)
{
  auto lhs = finite_vertices(before);
  auto rhs = finite_vertices(after);
  assert(lhs.size() == rhs.size());
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (!same_point(lhs[i], rhs[i])) return true;
  }
  return false;
}

template <
  typename C3t3T,
  typename ConcurrencyTag = CGAL::Mesh_smoothing_3::Parallel_if_available_tag
>
CGAL::Mesh_smoothing_3::Smoothing_status run_example_smoothing(C3t3T& c3t3, bool verbose, ConcurrencyTag tag = ConcurrencyTag{})
{
  auto status = CGAL::boundary_aware_mesh_smoothing<
    C3t3T,
    CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3T>
  >(
      c3t3,
      CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3T>(c3t3),
      CGAL::parameters::verbose(verbose).number_of_iterations(max_iterations).concurrency_tag(tag));
  assert_all_finite(c3t3);
  return status;
}

void assert_time_breakdown_close(CGAL::Mesh_smoothing_3::Smoothing_status const& status)
{
  const double sum = status.pre_processing_time + status.optimization_time;
  const double tol = (std::max)(1e-9, 1e-6 * (std::max)(1.0, std::abs(status.total_time)));
  assert(std::abs(status.total_time - sum) <= tol);
}

void test_verbose_has_no_impact()
{
  auto quiet = load_example_mesh();
  auto verbose = quiet;

  auto quiet_status = run_example_smoothing(quiet, false);
  auto verbose_status = run_example_smoothing(verbose, true);

  assert_structure_counts_preserved(quiet, verbose);
  assert(quiet_status.return_code == verbose_status.return_code);
  assert(quiet_status.nb_iterations == verbose_status.nb_iterations);
  assert(quiet_status.nb_vertex_updates == verbose_status.nb_vertex_updates);
  assert(quiet_status.nb_metric_evaluations == verbose_status.nb_metric_evaluations);
  assert(quiet_status.nb_initial_invalid_elements == verbose_status.nb_initial_invalid_elements);
  assert(quiet_status.nb_invalid_elements == verbose_status.nb_invalid_elements);
  assert_time_breakdown_close(quiet_status);
  assert_time_breakdown_close(verbose_status);
  assert_same_vertices(quiet, verbose);
}

void test_example_integration()
{
  auto before = load_example_mesh();
  auto after = before;

  auto status = run_example_smoothing(after, false);

  assert_structure_counts_preserved(before, after);
  assert(status.nb_initial_invalid_elements == 0);
  assert(status.nb_invalid_elements == 0);
  assert(status.valid_mesh());
  assert(status.nb_iterations > 0);
  assert(status.nb_vertex_updates > 0);
  assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::CONVERGENCE_REACHED);
  assert_time_breakdown_close(status);
  assert(any_vertex_moved(before, after));
}

void test_no_early_stopping()
{
  auto c3t3 = load_example_mesh();

  auto status = CGAL::boundary_aware_mesh_smoothing(
      c3t3,
      CGAL::Mesh_smoothing_3::C3t3_mesh_projector(c3t3),
      CGAL::parameters::verbose(false).number_of_iterations(2));

  assert(status.valid_mesh());
  assert(status.nb_iterations == 2);
  assert(status.nb_vertex_updates > 100);
  assert(status.return_code == CGAL::Mesh_smoothing_3::Smoothing_return_code::CONVERGENCE_REACHED); // sphere will converge in 2 iterations (to change)
  assert_time_breakdown_close(status);
}

void test_mono_core_reproducibility()
{

  auto first = load_example_mesh();
  auto second = load_example_mesh();

  auto first_status = run_example_smoothing<C3t3, CGAL::Sequential_tag>(first, false);
  auto second_status = run_example_smoothing<C3t3, CGAL::Sequential_tag>(second, false);

  assert_structure_counts_preserved(first, second);
  assert(first_status.return_code == second_status.return_code);
  assert(first_status.nb_iterations == second_status.nb_iterations);
  assert(first_status.nb_vertex_updates == second_status.nb_vertex_updates);
  assert(first_status.nb_metric_evaluations == second_status.nb_metric_evaluations);
  assert(first_status.nb_initial_invalid_elements == second_status.nb_initial_invalid_elements);
  assert(first_status.nb_invalid_elements == second_status.nb_invalid_elements);
  assert_time_breakdown_close(first_status);
  assert_time_breakdown_close(second_status);
  assert_same_vertices(first, second);
}

void test_random_inner_untangling_integration()
{
  std::mt19937 rng(424242u);
  auto target = load_example_mesh();
  auto projector = CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3>(target);
  auto perturbed = target;

  std::vector<Vertex_handle> interior;
  for (auto v : perturbed.triangulation().finite_vertex_handles()) {
    if (perturbed.in_dimension(v) == 3) interior.push_back(v);
  }
  assert(interior.size() >= 100);
  std::shuffle(interior.begin(), interior.end(), rng);

  CGAL::Bbox_3 bbox = perturbed.bbox();
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();
  const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double shift = untangle_amplitude * diag;
  std::uniform_real_distribution<double> du(-shift, shift);



  for (std::size_t i = 0; i < 100; ++i) {
    auto p = interior[i]->point();
    interior[i]->set_point(Point_3(CGAL::to_double(p.x()) + du(rng),
                                   CGAL::to_double(p.y()) + du(rng),
                                   CGAL::to_double(p.z()) + du(rng)));
  }


  auto status = CGAL::boundary_aware_mesh_smoothing(
      perturbed,
      projector,
      CGAL::parameters::verbose(false).number_of_iterations(untangle_iterations));


  assert(status.nb_initial_invalid_elements > 0);
  assert(status.nb_invalid_elements == 0);
  assert(status.valid_mesh());
  assert_time_breakdown_close(status);
}

void test_constraint_maps_integration()
{
  auto target = load_example_mesh();
  auto smoothed = target;
  auto& tr = smoothed.triangulation();

  // Find two different unmarked edges.
  Edge edge_forward;
  Edge edge_reverse;
  bool found_forward = false;
  bool found_reverse = false;

  auto edge_vertices = [](Edge const& e) {
    return std::make_pair(e.first->vertex(e.second), e.first->vertex(e.third));
  };

  for (auto const& e : tr.finite_edges()) {
    if (smoothed.is_in_complex(e))
      continue;

    if (!found_forward) {
      edge_forward = e;
      found_forward = true;
      continue;
    }

    const auto ev0 = edge_vertices(edge_forward);
    const auto ev1 = edge_vertices(e);

    // Prefer two distinct edges so both conventions are tested independently.
    if (ev1 != ev0 &&
        ev1 != std::make_pair(ev0.second, ev0.first)) {
      edge_reverse = e;
      found_reverse = true;
      break;
    }
  }

  assert(found_forward);
  assert(found_reverse);

  // Find an unmarked finite facet.
  Facet constrained_facet;
  bool found_facet = false;

  for (auto const& f : tr.finite_facets()) {
    if (!smoothed.is_in_complex(f)) {
      constrained_facet = f;
      found_facet = true;
      break;
    }
  }

  assert(found_facet);

  const auto e0 = edge_vertices(edge_forward);
  const auto e1 = edge_vertices(edge_reverse);
  const std::array<Vertex_handle, 3> fv = {
    constrained_facet.first->vertex((constrained_facet.second + 1) % 4),
    constrained_facet.first->vertex((constrained_facet.second + 2) % 4),
    constrained_facet.first->vertex((constrained_facet.second + 3) % 4)
  };

  std::map<std::pair<Vertex_handle, Vertex_handle>, bool> emap;

  // First edge in the iterator's orientation.
  emap[{e0.first, e0.second}] = true;

  // Second edge deliberately in reverse orientation.
  emap[{e1.second, e1.first}] = true;

  std::map<Facet, bool> fmap;

  // Also exercise the mirror representation.
  fmap[tr.mirror_facet(constrained_facet)] = true;

  std::vector<Vertex_handle> locked = {
    e0.first, e0.second,
    e1.first, e1.second,
    fv[0], fv[1], fv[2]
  };

  std::sort(locked.begin(), locked.end());
  locked.erase(std::unique(locked.begin(), locked.end()), locked.end());

  std::map<Vertex_handle, Point_3> before;
  for (auto v : tr.finite_vertex_handles())
    before[v] = v->point();

  auto projector =
    CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3>(target);

  auto status = CGAL::boundary_aware_mesh_smoothing(
    smoothed,
    projector,
    CGAL::parameters::
      edge_is_constrained_map(boost::make_assoc_property_map(emap))
      .facet_is_constrained_map(boost::make_assoc_property_map(fmap))
      .number_of_iterations(max_iterations));

  assert(status.nb_vertex_updates > 0);

  // Every vertex incident to a constrained edge/facet must be unchanged.
  for (auto v : locked) {
    assert(same_point(before[v], v->point()));
  }

  // Ensure the test did not pass simply because the entire mesh stayed fixed.
  bool unlocked_vertex_moved = false;

  for (auto v : tr.finite_vertex_handles()) {
    if (std::find(locked.begin(), locked.end(), v) != locked.end())
      continue;

    if (!same_point(before[v], v->point())) {
      unlocked_vertex_moved = true;
      break;
    }
  }

  assert(unlocked_vertex_moved);
  assert_all_finite(smoothed);
}

} // namespace

int main()
{
  test_verbose_has_no_impact();
  test_example_integration();
  test_no_early_stopping();
  test_mono_core_reproducibility();
  test_random_inner_untangling_integration();
  test_constraint_maps_integration();
  return EXIT_SUCCESS;
}
