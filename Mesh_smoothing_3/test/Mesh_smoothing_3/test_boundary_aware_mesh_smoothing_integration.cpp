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
  std::cerr << "loading_path=" << MS3_EXAMPLE_MESH_FILE << " open=" << is.good() << "\n";
  assert(is);
  const bool read_ok = CGAL::IO::read_MEDIT(is, c3t3.triangulation());
  std::cerr << "read_ok=" << read_ok
            << " vertices_after_read=" << c3t3.triangulation().number_of_vertices()
            << " cells_after_read=" << c3t3.triangulation().number_of_finite_cells() << "\n";
  assert(read_ok);
  c3t3.rescan_after_load_of_triangulation();
  std::cerr << "vertices_after_rescan=" << c3t3.triangulation().number_of_vertices()
            << " cells_after_rescan=" << c3t3.triangulation().number_of_finite_cells() << "\n";
  return c3t3;
}

void save_mesh(char const* filename, C3t3 const& c3t3)
{
  std::ofstream os(filename, std::ios::out);
  assert(os);
  CGAL::IO::write_MEDIT(os, c3t3.triangulation(), CGAL::parameters::all_vertices(true));
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

template <typename C3t3T>
CGAL::Mesh_smoothing_3::Smoothing_status run_example_smoothing(C3t3T& c3t3, bool verbose)
{
  auto status = CGAL::boundary_aware_mesh_smoothing(
      c3t3,
      CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3T>(c3t3),
      CGAL::parameters::verbose(verbose).number_of_iterations(max_iterations));
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
  std::cerr << "[test_verbose_has_no_impact]\n";
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
  std::cerr << "[test_example_integration]\n";
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

void test_mono_core_reproducibility()
{
  std::cerr << "[test_mono_core_reproducibility]\n";
  const int previous_threads = Eigen::nbThreads();
  Eigen::setNbThreads(1);
  assert(Eigen::nbThreads() == 1);

  auto first = load_example_mesh();
  auto second = load_example_mesh();

  auto first_status = run_example_smoothing(first, false);
  auto second_status = run_example_smoothing(second, false);

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

  Eigen::setNbThreads(previous_threads);
}

void test_random_inner_untangling_integration()
{
  std::cerr << "[test_random_inner_untangling_integration]\n";
  std::mt19937 rng(424242u);
  std::cerr << "mesh_path=" << MS3_EXAMPLE_MESH_FILE << '\n';
  auto target = load_example_mesh();
  std::cerr << "mesh_vertices=" << target.triangulation().number_of_vertices()
            << " cells=" << target.triangulation().number_of_finite_cells()
            << " facets=" << target.number_of_facets_in_complex()
            << " edges=" << target.number_of_edges_in_complex()
            << '\n';
  auto projector = CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3>(target);
  auto perturbed = target;
  save_mesh("test_boundary_aware_mesh_smoothing_integration_target.mesh", target);

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
  save_mesh("test_boundary_aware_mesh_smoothing_integration_perturbed.mesh", perturbed);

  auto status = CGAL::boundary_aware_mesh_smoothing(
      perturbed,
      projector,
      CGAL::parameters::verbose(true).number_of_iterations(untangle_iterations));
  save_mesh("test_boundary_aware_mesh_smoothing_integration_smoothed.mesh", perturbed);

  std::cerr << "initial_invalid=" << status.nb_initial_invalid_elements
            << " final_invalid=" << status.nb_invalid_elements
            << " valid=" << status.valid_mesh()
            << " it=" << status.nb_iterations
            << " updates=" << status.nb_vertex_updates
            << " return=" << int(status.return_code) << '\n';

  assert(status.nb_initial_invalid_elements > 0);
  assert(status.nb_invalid_elements == 0);
  assert(status.valid_mesh());
  assert_time_breakdown_close(status);
}

} // namespace

int main()
{
  test_verbose_has_no_impact();
  test_example_integration();
  test_mono_core_reproducibility();
  test_random_inner_untangling_integration();
  return EXIT_SUCCESS;
}
