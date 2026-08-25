// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include "test_common.h"

#include <cassert>
#include <cmath>

using namespace ms3_test;

namespace {

using Projection_mode = CGAL::Mesh_smoothing_3::Projection_weight_mode;

struct Projection_fixture
{
  C3t3 c3t3;
  std::array<Vertex_handle, 4> v{};
  Facet patch;
};

Projection_fixture make_projection_fixture()
{
  Projection_fixture f;
  auto& tr = f.c3t3.triangulation();

  f.v[0] = tr.insert(Point_3(0.0, 0.0, 0.30));
  f.v[1] = tr.insert(Point_3(1.0, 0.0, 0.30));
  f.v[2] = tr.insert(Point_3(0.0, 1.0, 0.30));
  f.v[3] = tr.insert(Point_3(0.0, 0.0, 1.40));

  auto cell = tr.finite_cells_begin();
  f.c3t3.add_to_complex(cell, 3);
  f.patch = facet_opposite(cell, f.v[3]);
  f.c3t3.add_to_complex(f.patch, 11);
  return f;
}

void test_direct_projector()
{
  auto fixture = make_fixture(false, false);
  CGAL::Mesh_smoothing_3::C3t3_mesh_projector<C3t3> projector(fixture.c3t3);

  auto patch = projector.patch_face_projection_plane(
      std::make_pair(11u, fixture.shared_facet),
      {fixture.v[0]->point(), fixture.v[1]->point(), fixture.v[2]->point()});
  assert(std::abs(CGAL::to_double(patch.origin().z())) < 1e-12);
  assert(std::abs(CGAL::to_double(patch.vector().z())) > 0.5);

  auto curve = projector.curve_edge_projection_line(
      std::make_pair(21u, fixture.curve_edge_1),
      {fixture.v[1]->point(), fixture.v[3]->point()});
  auto edge = fixture.v[3]->point() - fixture.v[1]->point();
  auto cross = CGAL::cross_product(curve.vector(), edge);
  assert(cross.squared_length() < 1e-12);
  auto d = curve.origin() - fixture.v[1]->point();
  assert(CGAL::cross_product(d, edge).squared_length() < 1e-12);
}

template <typename CTS>
double run_projection_case(CTS const& cts)
{
  auto fixture = make_projection_fixture();
  std::map<Vertex_handle, bool> vmap;
  vmap[fixture.v[3]] = true;

  auto status = CGAL::boundary_aware_mesh_smoothing(
      fixture.c3t3, cts,
      CGAL::parameters::vertex_is_constrained_map(boost::make_assoc_property_map(vmap))
          .number_of_iterations(1));

  (void)status;
  double sum = 0.;
  for (int i = 0; i < 3; ++i) {
    sum += std::abs(CGAL::to_double(fixture.v[i]->point().z()) - 0.25);
  }
  return sum;
}

void test_projection_weights()
{
  using GT = Recording_cts<C3t3>::Geom_traits;
  using Vector_3 = typename GT::Vector_3;
  using Tangent_space = typename Recording_cts<C3t3>::Tangent_space;

  Recording_cts<C3t3> none_cts;
  none_cts.patch_results.emplace(11, Tangent_space{Point_3(0., 0., 0.25), Vector_3(0., 0., 1.),
                                                   Projection_mode::NONE, 1.});
  double none_dist = run_projection_case(none_cts);

  Recording_cts<C3t3> zero_cts;
  zero_cts.patch_results.emplace(11, Tangent_space{Point_3(0., 0., 0.25), Vector_3(0., 0., 1.),
                                                   Projection_mode::CUSTOM, 0.0});
  double zero_dist = run_projection_case(zero_cts);

  Recording_cts<C3t3> soft_cts;
  soft_cts.patch_results.emplace(11, Tangent_space{Point_3(0., 0., 0.25), Vector_3(0., 0., 1.),
                                                   Projection_mode::SOFT, 1.});
  double soft_dist = run_projection_case(soft_cts);

  Recording_cts<C3t3> strong_cts;
  strong_cts.patch_results.emplace(11, Tangent_space{Point_3(0., 0., 0.25), Vector_3(0., 0., 1.),
                                                     Projection_mode::STRONG, 1.});
  double strong_dist = run_projection_case(strong_cts);

  assert(std::abs(none_dist - zero_dist) < 1e-12);
  assert(strong_dist <= soft_dist);
  assert(soft_dist <= none_dist + 1e-12);
}

void test_no_projection_smoke()
{
  auto fixture = make_fixture(false, false);
  auto before = finite_vertices(fixture.c3t3);

  auto status = CGAL::boundary_aware_mesh_smoothing(
      fixture.c3t3,
      CGAL::Mesh_smoothing_3::C3t3_no_projection<C3t3>{},
      CGAL::parameters::number_of_iterations(1));

  assert_all_finite(fixture.c3t3);
  assert(status.nb_vertex_updates > 0);
  auto after = finite_vertices(fixture.c3t3);
  assert(before.size() == after.size());
  bool any_moved = false;
  for (std::size_t i = 0; i < before.size(); ++i) {
    if (!same_point(before[i], after[i])) any_moved = true;
  }
  assert(any_moved);
  assert(status.valid_mesh() == (status.nb_invalid_elements == 0));
}

} // namespace

int main()
{
  test_direct_projector();
  test_projection_weights();
  test_no_projection_smoke();
  return EXIT_SUCCESS;
}
