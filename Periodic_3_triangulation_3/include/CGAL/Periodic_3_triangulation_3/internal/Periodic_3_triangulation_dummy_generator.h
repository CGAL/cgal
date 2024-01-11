// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#include <CGAL/license/Periodic_3_triangulation_3.h>

#ifdef CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H

template < class GT, class TDS >
inline std::vector<typename Periodic_3_triangulation_3<GT,TDS>::Vertex_handle >
Periodic_3_triangulation_3<GT, TDS>::
insert_generic_dummy_points()
{
  // the "info" field is the index of the neighbor offset in `neigh_offsets`
  using DT3_Vb = CGAL::Triangulation_vertex_base_with_info_3<std::size_t, GT>;
  using DT3_Cb = CGAL::Delaunay_triangulation_cell_base_3<GT>;
  using DT3_TDS = CGAL::Triangulation_data_structure_3<DT3_Vb, DT3_Cb>;
  using DT3 = CGAL::Delaunay_triangulation_3<GT, DT3_TDS>;
  using DT3_VH = typename DT3::Vertex_handle;
  using DT3_CH = typename DT3::Cell_handle;

  // Compute the number of subdivisions in all directions ------------------------------------------

  const std::array<FT, 3> spans = { domain().xmax() - domain().xmin(),
                                    domain().ymax() - domain().ymin(),
                                    domain().zmax() - domain().zmin() };

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 1)
  std::cout << "Domain spans:\n"
            << domain().xmin() << " " << domain().xmax() << "\n"
            << domain().ymin() << " " << domain().ymax() << "\n"
            << domain().zmin() << " " << domain().zmax() << std::endl;
#endif

  std::size_t min_pos = std::min_element(std::cbegin(spans), std::cend(spans)) - std::cbegin(spans);
  std::size_t max_pos = std::max_element(std::cbegin(spans), std::cend(spans)) - std::cbegin(spans);
  std::size_t mid_pos;

  if(min_pos == max_pos) // cubic
  {
    min_pos = 0; // x
    mid_pos = 1; // y
    max_pos = 2; // z
  }
  else
  {
    mid_pos = (min_pos + 1) % 3;
    if(mid_pos == max_pos)
      mid_pos = (max_pos + 1) % 3;

    if(min_pos > mid_pos && spans[min_pos] == spans[mid_pos])
      std::swap(min_pos, mid_pos); // just for convenience
  }

  CGAL_assertion(min_pos < 3 && mid_pos < 3 && max_pos < 3);
  CGAL_assertion(min_pos != max_pos && mid_pos != min_pos && mid_pos != max_pos);

  std::array<int, 3> nums_steps;
  std::array<FT, 3> steps;

  // Min:
  // this doesn't work for P3M3 + sharp features due to weights creating a harder constraint of 1/64th domain_size²
// #define CGAL_P3T3_USE_EXPERIMENTAL_LARGE_STEP_IN_DUMMY_GENERATION
#ifdef CGAL_P3T3_USE_EXPERIMENTAL_LARGE_STEP_IN_DUMMY_GENERATION
  nums_steps[min_pos] = 3;
#else
  nums_steps[min_pos] = 6;
#endif
  steps[min_pos] = spans[min_pos] / nums_steps[min_pos];

  // Mid: do not use the min step, but redistribute the error between nums_steps[mid_pos] * min_step and the actual span
  nums_steps[mid_pos] = int(to_interval(spans[mid_pos] / steps[min_pos]).first); // flooring, "min" is not a typo
  steps[mid_pos] = spans[mid_pos] / nums_steps[mid_pos];

  // Max: smaller step in the max span direction as to avoid cospherical configurations
#ifdef CGAL_P3T3_USE_EXPERIMENTAL_LARGE_STEP_IN_DUMMY_GENERATION
  const FT minor_step = spans[min_pos] / FT(4); // a ratio of 3:4 makes for nicely shaped tetrahedra
#else
  const FT minor_step = spans[min_pos] / FT(8); // a ratio of 6:8 makes for nicely shaped tetrahedra
#endif
  nums_steps[max_pos] = int(to_interval(spans[max_pos] / minor_step).first); // flooring
  CGAL_assertion(nums_steps[max_pos] >= steps[min_pos]);

  // Important! Consecutive levels in the max length have a shift (that's the `k % 2 != 0` part).
  // Hence, to get proper periodicity, the number of steps needs to be even, otherwise the level
  // `max - 1` and `0` would have the same shift and the Delaunay stars would not be identical everywhere.
  if(nums_steps[max_pos] % 2 != 0)
    ++nums_steps[max_pos];

  // now that the number of steps is known, re-adjust the step length so that it matches the height
  steps[max_pos] = spans[max_pos] / nums_steps[max_pos];

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 1)
  std::cout << "min|mid|max: " << min_pos << " " << mid_pos << " " << max_pos << std::endl;
  std::cout << "nums_steps[min_pos] = " << nums_steps[min_pos] << std::endl;
  std::cout << "steps[min_pos] = " << steps[min_pos] << std::endl;
  std::cout << "nums_steps[mid_pos] = " << nums_steps[mid_pos] << std::endl;
  std::cout << "steps[mid_pos] = " << steps[mid_pos] << std::endl;
  std::cout << "minor step = " << minor_step << std::endl;
  std::cout << "nums_steps[max_pos] = " << nums_steps[max_pos] << std::endl;
  std::cout << "steps[max_pos] = " << steps[max_pos] << std::endl;
#endif

  auto grid_offset_to_lattice_offset = [&](const Offset& loff) -> Offset
  {
    const int hloz = int(std::ceil(double(loff[max_pos]) / 2.));

    Offset goff;
    goff[min_pos] = loff[min_pos] + hloz;
    goff[mid_pos] = loff[mid_pos] + hloz;
    goff[max_pos] = loff[max_pos];
    return goff;
  };

  auto lattice_offset_to_grid_offset = [&](const Offset& loff) -> Offset
  {
    const int hloz = int(std::ceil(double(loff[max_pos]) / 2.));

    Offset goff;
    goff[min_pos] = loff[min_pos] - hloz;
    goff[mid_pos] = loff[mid_pos] - hloz;
    goff[max_pos] = loff[max_pos];
    return goff;
  };

  // This should be prefered from building points from the lattice offset
  // because the grid offset (by construction) aligns with the canonical domain,
  // and there is thus fewer numerical errors.
  auto construct_point_from_grid_offset = [&](const Offset& goff) -> Point_3
  {
    std::array<FT, 3> coords;
    coords[min_pos] = goff[min_pos] * steps[min_pos];
    coords[mid_pos] = goff[mid_pos] * steps[mid_pos];
    coords[max_pos] = goff[max_pos] * steps[max_pos];

    if(goff[max_pos] % 2 != 0)
    {
      coords[min_pos] += 0.5 * steps[min_pos];
      coords[mid_pos] += 0.5 * steps[mid_pos];
    }

    return { coords[0], coords[1], coords[2] };
  };

  auto construct_point_from_lattice_offset = [&](const Offset& loff) -> Point_3
  {
    // @fixme? verify that the lattice is indeed already reduced

    // For `min_pos = 0`, `mid_pos = 1`, `max_pos = 2`, this returns:
    //   CGAL::Origin + loff[0]*lv0 + loff[1]*lv1 + loff[2]*lv2
    // with
    //   Vector lv0 {        x_step,             0,      0 };
    //   Vector lv1 {             0,        y_step,      0 };
    //   Vector lv2 { -0.5 * x_step, -0.5 * y_step, z_step };

    std::array<FT, 3> coords;
    coords[min_pos] = loff[min_pos] * steps[min_pos] - loff[max_pos] * steps[min_pos] / FT(2);
    coords[mid_pos] = loff[mid_pos] * steps[mid_pos] - loff[max_pos] * steps[mid_pos] / FT(2);
    coords[max_pos] = loff[max_pos] * steps[max_pos];

    return { coords[0], coords[1], coords[2] };
  };

  CGAL_USE(construct_point_from_lattice_offset);

  // Create the dummy points -----------------------------------------------------------------------

  std::vector<Point_3> dummy_points;
  dummy_points.reserve(nums_steps[0] * nums_steps[1] * nums_steps[2]);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 32)
  DT3 dt3;
#endif

  std::array<FT, 3> coords;
  for(int i=0; i<nums_steps[min_pos]; ++i)
  {
    for(int j=0; j<nums_steps[mid_pos]; ++j)
    {
      for (int k=0; k<nums_steps[max_pos]; ++k)
      {
        coords[min_pos] = domain().min_coord(static_cast<int>(min_pos)) + i * steps[min_pos];
        coords[mid_pos] = domain().min_coord(static_cast<int>(mid_pos)) + j * steps[mid_pos];

        if(k % 2 != 0)
        {
          coords[min_pos] += steps[min_pos] / FT(2);
          coords[mid_pos] += steps[mid_pos] / FT(2);
        }

        coords[max_pos] = domain().min_coord(static_cast<int>(max_pos)) + k * steps[max_pos];

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 2)
        std::cout << "Add dummy: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
#endif
        dummy_points.emplace_back(coords[0], coords[1], coords[2]);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 32)
        DT3_VH dt3_v = dt3.insert(dummy_points.back());
#endif
      }
    }
  }

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 1)
  std::cout << dummy_points.size() << " dummy points" << std::endl;
#endif

  // Build the periodic triangulation --------------------------------------------------------------

  clear();
  tds().set_dimension(3);
  set_cover(CGAL::make_array(1,1,1));

  const std::size_t nv = static_cast<std::size_t>(nums_steps[0] * nums_steps[1] * nums_steps[2]);

  // Create vertices
  std::vector<Vertex_handle> vertices;
  vertices.reserve(nv);

  std::size_t id = 0;
  for(int i=0; i<nums_steps[min_pos]; ++i)
  {
    for(int j=0; j<nums_steps[mid_pos]; ++j)
    {
      for(int k=0; k<nums_steps[max_pos]; ++k)
      {
        Vertex_handle vh = tds().create_vertex();
        vertices.push_back(vh);
        vh->set_point(Point(dummy_points[id++])); // wrap with Point for regular triangulations
      }
    }
  }

  // Create cells
  std::vector<Offset> star_lnoffs;
  std::vector<std::array<std::size_t, 4> > star_cells;

  auto compute_periodic_star = [&]() -> void
  {
    // Lattice offsets sufficient to get the complete periodic star
    //
    // Considering all offsets is likely an overkill because it's a special lattice with bounded angles
    // and ratio, but it's a tiny triangulation so it does not matter
    const std::array<Offset, 75> neigh_offsets =
    {{
      {0, 0, 0},

      // 14 offsets that are entirely contained within the scaled domain
      {-1, -1, -1}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, -1, -1},
      {0, 1, 1}, {-1, 0, -1}, {-1, -1, 0}, {1, 1, 1}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0},

      // 36 offsets that have a guaranteed intersection with the scaled domain
      {1, 2, 0}, {1, 0, 2}, {-1, -2, -2}, {2, 1, 0}, {2, 0, 1}, {1, -1, -1}, {0, 1, 2},
      {-2, -1, -2}, {0, 2, 1}, {-1, 1, -1}, {-2, -2, -1}, {-1, -1, 1}, {1, 1, 2},
      {-1, -1, -2}, {1, 2, 1}, {0, 1, -1}, {-1, -2, -1}, {0, -1, 1}, {2, 1, 1},
      {1, 0, -1}, {1, -1, 0}, {-2, -1, -1}, {-1, 0, 1}, {-1, 1, 0}, {1, 2, 2},
      {-1, 0, -2}, {-1, -2, 0}, {2, 1, 2}, {0, -1, -2}, {2, 2, 1}, {1, 1, -1},
      {0, -2, -1}, {1, -1, 1}, {-2, -1, 0}, {-2, 0, -1}, {-1, 1, 1},

      // 24 offsets that might have an intersection with the scaled domain (6 of them will)
      {3, 2, 1}, {2, 1, -1}, {3, 1, 2}, {2, -1, 1}, {1, -1, -2}, {1, -2, -1},
      {2, 3, 1}, {1, 2, -1}, {1, 3, 2}, {-1, 2, 1}, {-1, 1, -2}, {-2, 1, -1},
      {2, 1, 3}, {1, -1, 2}, {1, 2, 3}, {-1, 1, 2}, {-1, -2, 1}, {-2, -1, 1},
      {-1, -2, -3}, {-1, -3, -2}, {-2, -1, -3}, {-3, -1, -2}, {-2, -3, -1}, {-3, -2, -1}
    }};

    star_lnoffs.reserve(neigh_offsets.size());
    star_lnoffs.push_back(neigh_offsets[0]);

    DT3 dt3;
    DT3_VH v0 = dt3.insert(CGAL::ORIGIN);
    CGAL_assertion(v0 != DT3_VH());
    v0->info() = 0;

    for(std::size_t oi=1; oi<neigh_offsets.size(); ++oi)
    {
      const Offset& lnoff = neigh_offsets[oi];

      Offset lidx;
      lidx[min_pos] = lnoff[0]; // not a typo, this aligns the neigh offsets with the min/mid/max
      lidx[mid_pos] = lnoff[1];
      lidx[max_pos] = lnoff[2];

      // Using gp because it's less likely to produce numerical errors
      const Offset gidx = lattice_offset_to_grid_offset(lidx);
      const Point_3 gp = construct_point_from_grid_offset(gidx);

      CGAL_warning_code(const Point_3 lp = construct_point_from_lattice_offset(lidx);)
      CGAL_warning_code(if(CGAL::squared_distance(gp, lp) > 1e-10 * spans[min_pos])) // might fail with inexact constructions
      CGAL_warning_code({ std::cout << "WARNING: gp/lp " << gp << " ||| " << lp << std::endl; })

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
      std::cout << "Lattice Offset: " << lidx[min_pos] << " " << lidx[mid_pos] << " " << lidx[max_pos]
                << " yields p: " << gp << std::endl;
#endif

      DT3_VH vh = dt3.insert(gp);
      CGAL_assertion(vh != DT3_VH() && vh != v0);
      vh->info() = oi;

      star_lnoffs.push_back(lidx);
    }

    std::list<DT3_CH> cells;
    dt3.incident_cells(v0, std::back_inserter(cells));

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
    std::set<Offset> neighbors_in_use;
#endif

    star_cells.reserve(cells.size());
    for(DT3_CH ch : cells)
    {
      // To avoid keeping tracks of cells being already created or not, we only insert cells
      // when they appear in the incident cells of the vertex with the smallest global ID within the call.
      //
      // Put the center vertex at pos 0 in the cell so that it is simple to find it
      constexpr std::array<std::array<int, 4>, 4> permutations =
      {{
        { 0, 1, 2, 3 },
        { 1, 2, 0, 3 },
        { 2, 0, 1, 3 },
        { 3, 1, 0, 2 },
      }};

      const int v0_pos = ch->index(v0);
      CGAL_postcondition(v0_pos >= 0 && v0_pos <= 3);
      CGAL_assertion(ch->vertex(permutations[v0_pos][0]) == v0);
      CGAL_assertion(ch->vertex(permutations[v0_pos][0])->info() == 0);

      star_cells.emplace_back(std::array<std::size_t, 4>{ ch->vertex(permutations[v0_pos][0])->info(),
                                                          ch->vertex(permutations[v0_pos][1])->info(),
                                                          ch->vertex(permutations[v0_pos][2])->info(),
                                                          ch->vertex(permutations[v0_pos][3])->info() });

      CGAL_assertion_code
      (
        for(int pos1=0; pos1<4; ++pos1)
          for(int pos2=0; pos2<4; ++pos2)
            if(pos1 != pos2)
              CGAL_assertion(star_cells.back()[pos1] != star_cells.back()[pos2]);
      )

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
      for(int i=0; i<4; ++i)
        neighbors_in_use.insert(neigh_offsets[ch->vertex(i)->info()]);
#endif
    }

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
    std::cout << neighbors_in_use.size() << " lattice offsets in the canonical star:\n";
    for(const Offset& loff : neighbors_in_use)
      std::cout << loff[min_pos] << " " << loff[mid_pos] << " " << loff[max_pos] << std::endl;

    std::cout << "as grid offsets:\n";
    for(const Offset& loff : neighbors_in_use)
    {
      const Offset goff = lattice_offset_to_grid_offset(loff);
      std::cout << goff[min_pos] << " " << goff[mid_pos] << " " << goff[max_pos] << std::endl;
    }
#endif
  };

  compute_periodic_star();

  // @todo might be able to avoid a map by numbering uniquely cells
  // with a clear vertex_id<->cell_id formula
  std::map<std::set<std::size_t>, std::pair<Facet, Facet> > neighboring_facets;

  for(int i=0; i<nums_steps[min_pos]; ++i)
  {
    for(int j=0; j<nums_steps[mid_pos]; ++j)
    {
      for(int k=0; k<nums_steps[max_pos]; ++k)
      {
        auto idx_to_global = [&](const Offset& v_idx) -> std::size_t
        {
          CGAL_assertion(v_idx[min_pos] >= 0 && v_idx[min_pos] < nums_steps[min_pos]);
          CGAL_assertion(v_idx[mid_pos] >= 0 && v_idx[mid_pos] < nums_steps[mid_pos]);
          CGAL_assertion(v_idx[max_pos] >= 0 && v_idx[max_pos] < nums_steps[max_pos]);

          return v_idx[min_pos] * nums_steps[mid_pos] * nums_steps[max_pos]
               + v_idx[mid_pos] * nums_steps[max_pos]
               + v_idx[max_pos];
        };

        auto add_facet = [&](const std::set<std::size_t>& face_vertices_ids, const Facet& f) -> void
        {
          CGAL_assertion(f.first != Cell_handle() && f.second >= 0 && f.second <= 3);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
          std::cout << std::endl;
          std::cout << "Add facet";
          for(const std::size_t vi : face_vertices_ids)
            std::cout << " " << vi;
          std::cout << std::endl;

          std::cout << "Positions\n";
          for(const std::size_t vi : face_vertices_ids)
            std::cout << vertices[vi]->point() << std::endl;
          std::cout << std::endl;
#endif

          auto insertion_result = neighboring_facets.emplace(face_vertices_ids, std::make_pair(f, Facet()));
          if(!insertion_result.second) // already exists in the facet map
          {
            std::pair<Facet, Facet>& fp = insertion_result.first->second;
            CGAL_assertion(fp.first.first != Cell_handle()); // one existing neighbor since insertion failed
            CGAL_assertion(fp.second.first == Cell_handle()); // two neighbors at most
            CGAL_assertion(fp.second.first != f.first); // neighbors are different
            fp.second = f;

            CGAL_assertion(fp.second.first != Cell_handle()); // properly set up

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
            std::cout << "facets" << std::endl;
            std::cout << &*(fp.first.first->vertex((fp.first.second+1)%4)) << std::endl;
            std::cout << &*(fp.first.first->vertex((fp.first.second+2)%4)) << std::endl;
            std::cout << &*(fp.first.first->vertex((fp.first.second+3)%4)) << std::endl;
            std::cout << &*(fp.second.first->vertex((fp.second.second+1)%4)) << std::endl;
            std::cout << &*(fp.second.first->vertex((fp.second.second+2)%4)) << std::endl;
            std::cout << &*(fp.second.first->vertex((fp.second.second+3)%4)) << std::endl;
            std::cout << fp.first.first->vertex((fp.first.second+1)%4)->point() << std::endl;
            std::cout << fp.first.first->vertex((fp.first.second+2)%4)->point() << std::endl;
            std::cout << fp.first.first->vertex((fp.first.second+3)%4)->point() << std::endl;
            std::cout << fp.second.first->vertex((fp.second.second+1)%4)->point() << std::endl;
            std::cout << fp.second.first->vertex((fp.second.second+2)%4)->point() << std::endl;
            std::cout << fp.second.first->vertex((fp.second.second+3)%4)->point() << std::endl;

            std::set<Vertex_handle> f1_vs = { fp.first.first->vertex((fp.first.second+1)%4),
                                              fp.first.first->vertex((fp.first.second+2)%4),
                                              fp.first.first->vertex((fp.first.second+3)%4) };
            std::set<Vertex_handle> f2_vs = { fp.second.first->vertex((fp.second.second+1)%4),
                                              fp.second.first->vertex((fp.second.second+2)%4),
                                              fp.second.first->vertex((fp.second.second+3)%4) };
            CGAL_assertion(f1_vs == f2_vs);
#endif

            // neighbors could be set here instead of storing pairs and doing it later,
            // but it doesn't cost much and it enables checking that every face has exactly
            // two incident cells
          }

          const std::pair<Facet, Facet>& fp = insertion_result.first->second;
          CGAL_USE(fp);
          CGAL_assertion(fp.first.first != Cell_handle()); // properly set up
        };

        auto add_cell = [&](const std::array<std::size_t, 4>& star_cell) -> Cell_handle
        {
#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
          std::cout << std::endl;
          std::cout << "-- Considering cell " << star_cell[0] << " " << star_cell[1] << " " << star_cell[2] << " " << star_cell[3] << std::endl;
#endif

          // convert to lattice offset, add (lattice) neighbor offset, convert back to grid

          Offset idx;
          idx[min_pos] = i;
          idx[mid_pos] = j;
          idx[max_pos] = k;

          const Offset lidx = grid_offset_to_lattice_offset(idx);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
          std::cout << "(Recall that min/mid/max are: " << min_pos << " " << mid_pos << " " << max_pos << ")" << std::endl;
          std::cout << "(Recall that steps are: " << steps[min_pos] << " " << steps[mid_pos] << " " << steps[max_pos] << ")" << std::endl;
          std::cout << "(Recall that nums steps are: " << nums_steps[min_pos] << " " << nums_steps[mid_pos] << " " << nums_steps[max_pos] << ")" << std::endl;
          std::cout << "grid idx: " << idx[min_pos] << " " << idx[mid_pos] << " " << idx[max_pos] << std::endl;
          std::cout << "lattice idx: " << lidx[min_pos] << " " << lidx[mid_pos] << " " << lidx[max_pos] << std::endl;
#endif

          std::array<std::size_t, 4> glob_idx;
          std::array<CGAL::Periodic_3_offset_3, 4> offsets { };

          for(int pos=0; pos<4; ++pos)
          {
            const Offset& lnoff = star_lnoffs[star_cell[pos]];

            // lnoff has already been shuffled to match min/mid/max
            Offset lnidx = lidx + lnoff;

            // canonical (to be) ijk of the neighbor
            Offset cidx = lattice_offset_to_grid_offset(lnidx);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
            std::cout << std::endl;
            std::cout << "lattice neighbor offset[" << pos << "] = " << lnoff[min_pos] << " " << lnoff[mid_pos] << " " << lnoff[max_pos] << std::endl;
            std::cout << "lattice neighbor idx[" << pos << "] = " << lnidx[min_pos] << " " << lnidx[mid_pos] << " " << lnidx[max_pos] << std::endl;
            std::cout << "canonical idx[" << pos << "] = " << cidx[min_pos] << " " << cidx[mid_pos] << " " << cidx[max_pos] << " (before)" << std::endl;
#endif

            CGAL_warning_code(const Point_3 gp = construct_point_from_grid_offset(cidx);)
            CGAL_warning_code(const Point_3 lp = construct_point_from_lattice_offset(lnidx);)
            CGAL_warning_code(if(CGAL::squared_distance(gp, lp) > 1e-10 * spans[min_pos])) // might fail with inexact constructions
            CGAL_warning_code({ std::cout << "WARNING: gp/lp " << gp << " ||| " << lp << std::endl; })

            // canonicalize using indices to avoid constructions
            for(int l=0; l<3; ++l)
            {
              offsets[pos][l] = 0;

              while(cidx[l] >= nums_steps[l])
              {
                cidx[l] -= nums_steps[l];
                offsets[pos][l] += 1;
              }

              while(cidx[l] < 0)
              {
                cidx[l] += nums_steps[l];
                offsets[pos][l] -= 1;
              }
            }

            glob_idx[pos] = idx_to_global(cidx);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
            std::cout << "final canon idx[" << pos << "] = " << cidx[min_pos] << " " << cidx[mid_pos] << " " << cidx[max_pos] << std::endl;
            std::cout << "periodic offset[" << pos << "] = " << offsets[pos] << std::endl;
            std::cout << "Glob ID[" << pos << "] = " << glob_idx[pos]<< std::endl;
            std::cout << "Position[" << pos << "] = " << vertices[glob_idx[pos]]->point() << std::endl;
#endif

            CGAL_assertion(glob_idx[pos] < nv);
          }

          CGAL_assertion_code
          (
            for(int pos1=0; pos1<4; ++pos1)
              for(int pos2=0; pos2<4; ++pos2)
                if(pos1 != pos2)
                  CGAL_assertion(glob_idx[pos1] != glob_idx[pos2]);
          )

          // trick to add facets only once per cell; can't really do that for all operations
          // due to cells on the border (or, it could be done, but then
          // is_not_on_top_right_up_corner must be reworked)
          //
          // The trick works because glob_idx[0] is always the center vertex
          if(glob_idx[0] != (std::min)({glob_idx[0], glob_idx[1], glob_idx[2], glob_idx[3]}))
            return Cell_handle();

          // here and below, it's a new cell

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
          std::cout << "\nCell position around star center" << std::endl;
          std::cout << construct_point(vertices[glob_idx[0]]->point(), offsets[0]) << std::endl;
          std::cout << construct_point(vertices[glob_idx[1]]->point(), offsets[1]) << std::endl;
          std::cout << construct_point(vertices[glob_idx[2]]->point(), offsets[2]) << std::endl;
          std::cout << construct_point(vertices[glob_idx[3]]->point(), offsets[3]) << std::endl;
#endif

          CGAL_assertion(orientation(construct_point(vertices[glob_idx[0]]->point(), offsets[0]),
                                     construct_point(vertices[glob_idx[1]]->point(), offsets[1]),
                                     construct_point(vertices[glob_idx[2]]->point(), offsets[2]),
                                     construct_point(vertices[glob_idx[3]]->point(), offsets[3])) == CGAL::POSITIVE);

          Cell_handle ch = tds().create_cell();
          CGAL_assertion(ch != Cell_handle());

          CGAL_assertion_code(for(int pos=0; pos<4; ++pos) {)
          CGAL_assertion(vertices.at(glob_idx[pos]) != Vertex_handle());
          CGAL_assertion_code(})

          ch->set_vertices(vertices[glob_idx[0]], vertices[glob_idx[1]], vertices[glob_idx[2]], vertices[glob_idx[3]]);
          set_offsets(ch, offsets[0], offsets[1], offsets[2], offsets[3]);

          for(int pos=0; pos<4; ++pos)
            vertices[glob_idx[pos]]->set_cell(ch);

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 16)
          std::cout << "\nCell position in P3T3" << std::endl;
          std::cout << point(ch, 0) << std::endl;
          std::cout << point(ch, 1) << std::endl;
          std::cout << point(ch, 2) << std::endl;
          std::cout << point(ch, 3) << std::endl;
#endif

          // add faces to the neighboring map
          add_facet(std::set<std::size_t>{glob_idx[0], glob_idx[1], glob_idx[2]}, Facet(ch, 3));
          add_facet(std::set<std::size_t>{glob_idx[0], glob_idx[1], glob_idx[3]}, Facet(ch, 2));
          add_facet(std::set<std::size_t>{glob_idx[0], glob_idx[2], glob_idx[3]}, Facet(ch, 1));
          add_facet(std::set<std::size_t>{glob_idx[1], glob_idx[2], glob_idx[3]}, Facet(ch, 0));

          return ch;
        };

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 4)
        std::cout << std::endl;
        std::cout << " ====================== VERTEX " << i << " " << j << " " << k << " ==============" << std::endl;
#endif

        for(const auto& star_cell : star_cells)
          add_cell(star_cell);
      }
    }
  }

  // build neighboring info
  for(const auto& e : neighboring_facets)
  {
    const std::pair<Facet, Facet>& fp = e.second;
    const Cell_handle ch_a = fp.first.first, ch_b = fp.second.first;
    const int s_a = fp.first.second, s_b = fp.second.second;

    CGAL_assertion(ch_a != Cell_handle() && ch_b != Cell_handle() && ch_a != ch_b);
    CGAL_assertion(ch_a->neighbor(s_a) == Cell_handle() && ch_b->neighbor(s_b) == Cell_handle());

    ch_a->set_neighbor(s_a, ch_b);
    ch_b->set_neighbor(s_b, ch_a);
  }

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 8)
  std::cout << "facet map size: " << neighboring_facets.size() << std::endl;
  std::cout << number_of_vertices() << " nv" << std::endl;
  std::cout << number_of_cells() << " nc" << std::endl;
#endif

#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 2)
  CGAL_postcondition(is_valid(true));
#else
  CGAL_postcondition(is_valid());
#endif

  return vertices;
}

#endif // CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H
