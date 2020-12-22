// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_EXPERIMENTAL_H
#define CGAL_KSR_3_EXPERIMENTAL_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL {
namespace KSR_3 {

class Experimental {

  const std::vector<PVertex> merge_pvertices_on_ivertex(
    const FT min_time, const FT max_time,
    const PVertex& event_pvertex,
    const IVertex& ivertex,
    std::vector<PVertex>& pvertices,
    std::vector<IEdge>& crossed) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** merging " << str(event_pvertex) << " on " << str(ivertex) << std::endl;
    }

    crossed.clear();
    const KSR::size_t support_plane_idx = pvertices.front().first;
    const PVertex prev = pvertices.front();
    const PVertex next = pvertices.back();

    IEdge prev_iedge = null_iedge(), next_iedge = null_iedge();
    // std::cout << "starting from: " << segment_3(iedge(pvertices[1])) << std::endl;

    if (m_verbose) {
      std::cout << "- start from: " <<
      str(iedge(pvertices[1])) << " " << segment_3(iedge(pvertices[1])) << std::endl;
    }

    // Copy front/back to remember position/direction.
    PVertex front, back;
    if (pvertices.size() < 3) {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    } else if (pvertices.size() == 3 || pvertices.size() == 4) {

      // BUG: In this case, the point that is duplicated twice is not always copied.
      // To fix it, we copy the second point not from the original vertex but from the first
      // copy of that vertex.

      const auto& initial = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial.second));
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(front.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(front.second));

    } else if (pvertices.size() >= 5) {

      const auto& initial1 = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial1.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial1.second));

      const auto& initial2 = pvertices[pvertices.size() - 2];
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial2.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(initial2.second));

    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    }

    // auto pvertex_to_point =
    //   [&](const PVertex& a) -> Point_2 {
    //     return point_2(a);
    //   };

    // PFace fprev = pface_of_pvertex(prev);
    // Point_2 pprev = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fprev).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fprev).end(), pvertex_to_point));
    // PFace fnext = pface_of_pvertex(next);
    // Point_2 pnext = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fnext).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fnext).end(), pvertex_to_point));

    bool was_swapped = false;
    // if (CGAL::orientation(pprev, point_2(support_plane_idx, ivertex), pnext) == CGAL::LEFT_TURN) {
    //   std::cout << "swapped" << std::endl;
    //   was_swapped = true;
    //   std::swap(prev, next);
    //   std::swap(front, back);
    // }

    // Freeze vertices.
    for (std::size_t i = 1; i < pvertices.size() - 1; ++i) {
      PVertex& pvertex = pvertices[i];
      Point_2 point = point_2(support_plane_idx, ivertex);
      support_plane(pvertex).direction(pvertex.second) = CGAL::NULL_VECTOR;
      support_plane(pvertex).set_point(pvertex.second, point);
    }

    PVertex pvertex = pvertices[1];
    connect (pvertex, ivertex);

    if (m_verbose) {
      std::cout << "- frozen pvertex: " << str(pvertex) << std::endl;
    }
    // std::cout << point_3(pvertex) << std::endl;
    // std::cout << "removed pvertices:";

    // Join vertices.
    for (std::size_t i = 2; i < pvertices.size() - 1; ++ i)
    {
      // std::cout << " " << str(pvertices[i]) << std::endl;
      // std::cout << point_3(pvertices[i]) << std::endl;

      const auto he = mesh(support_plane_idx).halfedge(pvertices[i].second, pvertex.second);
      disconnect_ivertex (pvertices[i]);
      CGAL::Euler::join_vertex(he, mesh(support_plane_idx));
    }
    // std::cout << std::endl;

    auto i_iedges = incident_iedges(ivertex);
    std::vector<std::pair<IEdge, Direction_2> > iedges;
    std::copy (i_iedges.begin(), i_iedges.end(),
               boost::make_function_output_iterator
               ([&](const IEdge& ie) -> void
                {
                  if (intersected_planes(ie).find (support_plane_idx)
                      == intersected_planes(ie).end())
                    return;

                  Direction_2 dir (point_2 (support_plane_idx, opposite (ie, ivertex))
                                   - point_2 (support_plane_idx, ivertex));
                  iedges.push_back (std::make_pair (ie, dir));
                }));

    std::sort (iedges.begin(), iedges.end(),
               [&](const std::pair<IEdge, Direction_2>& a,
                   const std::pair<IEdge, Direction_2>& b) -> bool
               {
                 return a.second < b.second;
               });
    CGAL_assertion(iedges.size() != 0);

    if (m_verbose) {
      std::cout << "- neighbors: " << std::endl <<
      "prev = " << point_3(prev)  << " / " << direction(prev)  << std::endl <<
      "fron = " << point_3(front) << " / " << direction(front) << std::endl <<
      "back = " << point_3(back)  << " / " << direction(back)  << std::endl <<
      "next = " << point_3(next)  << " / " << direction(next)  << std::endl;
    }

    // std::cout << (iedge(next) != null_iedge()) << std::endl;
    // std::cout << "source: " << point_3(source(iedge(next))) << std::endl;
    // std::cout << "target: " << point_3(target(iedge(next))) << std::endl;
    // std::cout << "ivertex: " << point_3(ivertex) << std::endl;

    bool back_constrained = false;
    if ((iedge(next) != null_iedge()
         && (source(iedge(next)) == ivertex || target(iedge(next)) == ivertex))
        || (this->ivertex(next) != null_ivertex()
            && is_iedge (this->ivertex(next), ivertex)))
      back_constrained = true;

    bool front_constrained = false;
    if ((iedge(prev) != null_iedge()
         && (source(iedge(prev)) == ivertex || target(iedge(prev)) == ivertex))
        || (this->ivertex(prev) != null_ivertex()
            && is_iedge (this->ivertex(prev), ivertex)))
      front_constrained = true;

    std::vector<PVertex> new_vertices;
    if (back_constrained && front_constrained) // Closing case
    {
      if (m_verbose) {
        std::cout << "*** CLOSING CASE" << std::endl;
      }
    }
    else if (back_constrained) // Border case
    {
      if (m_verbose) {
        std::cout << "*** BACK BORDER CASE" << std::endl;
      }

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(prev) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const Direction_2 tmp_dir(tmp_prev - point_2(pvertex.first, ivertex));
      // std::cout << "tmp_dir: " << to_3d(prev.first, tmp_prev) << std::endl;

      std::reverse(iedges.begin(), iedges.end());

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (tmp_dir.counterclockwise_in_between(
          iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        const bool limit_reached = (line_idx(iedge) == other_side_limit);
        if (m_verbose) {
          std::cout << "- limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;
        }

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: BACK STD WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      std::vector<IEdge> all_crossed;
      // iedge_idx = first_idx; iter = 0;
      // while (true) {
      //   const IEdge& iedge = iedges[iedge_idx].first;
      //   bool limit_reached, bbox_reached;
      //   std::tie(limit_reached, bbox_reached) = is_occupied(pvertex, iedge);
      //   all_crossed.push_back(iedge);
      //   if (limit_reached || bbox_reached) {
      //     break;
      //   }
      //   iedge_idx = (iedge_idx + 1) % iedges.size();
      //   if (iter == 100) {
      //     CGAL_assertion_msg(false, "ERROR: BACK LIMIT WHY SO MANY ITERATIONS?");
      //   } ++iter;
      // }

      CGAL_assertion(crossed.size() != 0);
      if (m_verbose) {
        std::cout << "- crossed " << crossed.size() << " iedges:" << std::endl;
        for (const auto& iedge : crossed) {
          std::cout << segment_3(iedge) << std::endl;
        }
      }
      // CGAL_assertion(crossed[0] == all_crossed[0]);

      std::vector<Point_2> future_points;
      std::vector<Vector_2> future_directions;

      if (crossed.size() > all_crossed.size()) {
        future_points.resize(crossed.size());
        future_directions.resize(crossed.size());
        for (std::size_t i = 0; i < crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, back, prev, crossed[i], future_points[i], future_directions[i]);
          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
              prev_iedge = crossed[i];
            }
          }
        }
      } else {
        future_points.resize(all_crossed.size());
        future_directions.resize(all_crossed.size());
        for (std::size_t i = 0; i < all_crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, back, prev, all_crossed[i], future_points[i], future_directions[i]);
          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, prev, all_crossed[i])) {
              prev_iedge = all_crossed[i];
            }
          }
        }
      }

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "back saved: " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            i, back, prev, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          if (m_verbose) {
            std::cout << "- cropping" << std::endl;
          }
          PVertex cropped;
          if (prev_iedge != null_iedge() && prev_iedge == crossed[i]) {
            if (m_verbose) std::cout << "- prev parallel case" << std::endl;

            cropped = prev;
            Point_2 future_point; Vector_2 future_direction;
            const auto pair = this->border_prev_and_next(prev);
            const auto pprev = pair.first;
            compute_future_point_and_direction(
              i, prev, pprev, prev_iedge, future_point, future_direction);
            future_points[i] = future_point;
            future_directions[i] = future_direction;

          } else {
            if (m_verbose) std::cout << "- standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
            // future_point = future_points[i];
            // future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          // std::cout << "cropped point -> direction: " << point_2(cropped) << " -> " << direction(cropped) << std::endl;
          if (m_verbose) std::cout << "- cropped: " << point_3(cropped) << std::endl;
        }
        else // create triangle face
        {
          CGAL_assertion_msg(i == 1,
          "TODO: BACK, CAN WE HAVE MORE THAN 1 NEW FACE? IF YES, I SHOULD CHECK K FOR EACH!");
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);
          if (m_verbose) {
            std::cout << "- is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;
          }

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            if (m_verbose) std::cout << "- stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            if (m_verbose) std::cout << "- stop k" << std::endl;
            break;
          }

          // Create a new face.
          if (m_verbose) std::cout << "- adding new pface" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            if (m_verbose) std::cout << "- continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            if (m_verbose) std::cout << "- continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          if (m_verbose) {
            // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
            std::cout << "- k intersections after: " << this->k(pface) << std::endl;
          }

          const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          if (m_verbose) std::cout << "- propagated: " << point_3(propagated) << std::endl;
          const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }

      if (crossed.size() < all_crossed.size()) {
        const std::size_t csize = crossed.size();
        const std::size_t asize = all_crossed.size();

        if (m_verbose) {
          std::cout << "- crossed size: " << csize << std::endl;
          std::cout << "- all_crossed size: " << asize << std::endl;
        }

        const std::size_t num_extra_faces = asize - csize;
        CGAL_assertion(num_extra_faces > 0);
        if (num_extra_faces == 1) {

          if (m_verbose) std::cout << "- adding extra face" << std::endl;
          PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed.back());
          if (propagated == null_pvertex()) {
            CGAL_assertion_msg(false, "TODO: BACK, NULL PROPAGATED CASE!");
          } else {

            // std::cout << "propagated: " << point_3(propagated) << std::endl;
            CGAL_assertion(num_extra_faces == 1);

            // Old code.
            // const PFace pface = pface_of_pvertex(pvertex);
            // const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
            // this->k(new_pface) = this->k(pface);
            // previous = propagated;
            // continue;

            // New code!
            const auto opposite = this->opposite(all_crossed.back(), ivertex);
            const auto& mesh = this->mesh(pvertex);
            auto he = mesh.halfedge(pvertex.second, propagated.second);
            const PEdge qpedge(pvertex.first, mesh.edge(he));
            // std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

            PFace target_pface = null_pface();
            for (const auto pface : pfaces(pvertex.first)) {
              for (const auto pedge : pedges_of_pface(pface)) {
                if (pedge == qpedge) {
                  target_pface = pface;
                  break;
                }
              }
            }
            CGAL_assertion(target_pface != null_pface());

            const auto tt = pedges_of_pface(target_pface);
            std::vector<PEdge> pedges;
            pedges.reserve(tt.size());
            for (const auto t : tt) pedges.push_back(t);

            PEdge other_pedge = null_pedge();
            for (std::size_t j = 0; j < pedges.size(); ++j) {
              if (pedges[j] == qpedge) {
                const std::size_t jp = (j + 1) % pedges.size();
                const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
                const auto& pedge1 = pedges[jm];
                const auto& pedge2 = pedges[jp];
                const auto iv1 = this->ivertex(this->target(pedge1));
                const auto iv2 = this->ivertex(this->source(pedge2));
                if (iv1 == opposite) {
                  CGAL_assertion(iv2 != opposite);
                  other_pedge = pedge1;
                  break;
                } else if (iv2 == opposite) {
                  CGAL_assertion(iv1 != opposite);
                  other_pedge = pedge2;
                  break;
                } else {
                  CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
                }
              }
            }
            CGAL_assertion(other_pedge != null_pedge());
            // std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

            IEdge other_iedge;
            const auto& iedges = support_plane(pvertex).iedges();
            CGAL_assertion(has_iedge(other_pedge));
            const auto query_iedge = this->iedge(other_pedge);
            for (const auto& iedge : iedges) {
              if (iedge == query_iedge) continue;
              if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
                if (line_idx(query_iedge) == line_idx(iedge)) {
                  other_iedge = iedge;
                }
              }
            }
            CGAL_assertion(other_iedge != null_iedge());
            // std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

            CGAL_assertion(m_points.find(std::make_pair(pvertex.first, other_iedge)) != m_points.end());
            CGAL_assertion(m_directions.find(std::make_pair(pvertex.first, other_iedge)) != m_directions.end());
            const Point_2 future_point = m_points.at(std::make_pair(pvertex.first, other_iedge));
            const Vector_2 future_direction = m_directions.at(std::make_pair(pvertex.first, other_iedge));

            auto tmp = future_direction;
            tmp = KSR::normalize(tmp);
            // std::cout << "future tmp: " << to_3d(pvertex.first, point_2(propagated) + m_current_time * tmp) << std::endl;

            const auto before = propagated;
            propagated = add_pvertex(propagated.first, future_point);
            direction(propagated) = future_direction;
            new_vertices.push_back(propagated);

            // std::cout << "before: " << point_3(before) << std::endl;
            // std::cout << "propagated: " << point_3(propagated) << std::endl;

            const PFace pface = pface_of_pvertex(pvertex);
            const PFace new_pface = add_pface(std::array<PVertex, 4>{pvertex, before, propagated, previous});
            this->k(new_pface) = this->k(pface);
            previous = propagated;

            // CGAL_assertion_msg(false, "DEBUG THIS CASE!");

            const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(before.second, propagated.second));
            connect(pedge, other_iedge);
            connect(propagated, other_iedge);
            crossed.push_back(other_iedge);
          }
        } else {
          CGAL_assertion_msg(false, "TODO: BACK, CROSSED < LIMIT, MULTIPLE FACES!");
        }
      }

      if (crossed.size() == all_crossed.size()) {
        // continue...
      }

      if (crossed.size() > all_crossed.size()) {
        // continue ..
        // std::cout << "crossed size: " << crossed.size() << std::endl;
        // std::cout << "all crossed size: " << all_crossed.size() << std::endl;
        CGAL_assertion_msg(false, "TODO: BACK CROSSED > LIMIT!");
      }
    }
    else if (front_constrained) // Border case
    {
      if (m_verbose) {
        std::cout << "*** FRONT BORDER CASE" << std::endl;
      }

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(next) - point_2(pvertex));

      const FT next_time = last_event_time(next);
      // std::cout << next_time << std::endl;
      // std::cout << m_current_time << std::endl;
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(next_time >= FT(0));

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 tmp_dir(tmp_next - point_2(pvertex.first, ivertex));
      // std::cout << "tmp_dir: " << to_3d(next.first, tmp_next) << std::endl;

      if (was_swapped) {
        std::reverse(iedges.begin(), iedges.end());
      }

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (!was_swapped) {
          if (tmp_dir.counterclockwise_in_between(
            iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        } else {
          if (tmp_dir.counterclockwise_in_between(
            iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        }
      }

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        const bool limit_reached = (line_idx(iedge) == other_side_limit);

        if (m_verbose) {
          std::cout << "- limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;
        }

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: FRONT STD WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      std::vector<IEdge> all_crossed;
      // iedge_idx = first_idx; iter = 0;
      // while (true) {
      //   const IEdge& iedge = iedges[iedge_idx].first;
      //   bool limit_reached, bbox_reached;
      //   std::tie(limit_reached, bbox_reached) = is_occupied(pvertex, iedge);
      //   all_crossed.push_back(iedge);
      //   if (limit_reached || bbox_reached) {
      //     break;
      //   }
      //   iedge_idx = (iedge_idx + 1) % iedges.size();
      //   if (iter == 100) {
      //     CGAL_assertion_msg(false, "ERROR: FRONT LIMIT WHY SO MANY ITERATIONS?");
      //   } ++iter;
      // }

      CGAL_assertion(crossed.size() != 0);
      if (m_verbose) {
        std::cout << "- crossed " << crossed.size() << " iedges:" << std::endl;
        for (const auto& iedge : crossed) {
          std::cout << segment_3(iedge) << std::endl;
        }
      }
      // CGAL_assertion(crossed[0] == all_crossed[0]);

      std::vector<Point_2> future_points;
      std::vector<Vector_2> future_directions;

      if (crossed.size() > all_crossed.size()) {
        future_points.resize(crossed.size());
        future_directions.resize(crossed.size());
        for (std::size_t i = 0; i < crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, front, next, crossed[i], future_points[i], future_directions[i]);

          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
              next_iedge = crossed[i];
            }
          }
        }
      } else {
        future_points.resize(all_crossed.size());
        future_directions.resize(all_crossed.size());
        for (std::size_t i = 0; i < all_crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, front, next, all_crossed[i], future_points[i], future_directions[i]);

          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, next, all_crossed[i])) {
              next_iedge = all_crossed[i];
            }
          }
        }
      }

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "front saved: " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            i, front, next, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          if (m_verbose) std::cout << "- cropping" << std::endl;
          PVertex cropped;
          if (next_iedge != null_iedge() && next_iedge == crossed[i]) {
            if (m_verbose) std::cout << "- next parallel case" << std::endl;

            cropped = next;
            Point_2 future_point; Vector_2 future_direction;
            const auto pair = this->border_prev_and_next(next);
            const auto nnext = pair.second;
            compute_future_point_and_direction(
              i, next, nnext, next_iedge, future_point, future_direction);
            future_points[i] = future_point;
            future_directions[i] = future_direction;

          } else {
            if (m_verbose) std::cout << "- standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
            // future_point = future_points[i];
            // future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          CGAL_assertion(cropped != pvertex);
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          // std::cout << "cropped point -> direction: " << point_2(cropped) << " -> " << direction(cropped) << std::endl;
          if (m_verbose) std::cout << "- cropped: " << point_3(cropped) << std::endl;
        }
        else // create triangle face
        {
          CGAL_assertion_msg(i == 1,
          "TODO: FRONT, CAN WE HAVE MORE THAN 1 NEW FACE? IF YES, I SHOULD CHECK K FOR EACH!");
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);

          if (m_verbose) {
            std::cout << "- is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;
          }

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            if (m_verbose) std::cout << "- stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            if (m_verbose) std::cout << "- stop k" << std::endl;
            break;
          }

          // Create a new face.
          if (m_verbose) std::cout << "- adding new pface" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            if (m_verbose) std::cout << "- continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            if (m_verbose) std::cout << "- continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          if (m_verbose) {
            // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
            std::cout << "- k intersections after: " << this->k(pface) << std::endl;
          }

          const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          if (m_verbose) std::cout << "- propagated: " << point_3(propagated) << std::endl;
          const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }

      if (crossed.size() < all_crossed.size()) {
        const std::size_t csize = crossed.size();
        const std::size_t asize = all_crossed.size();

        if (m_verbose) {
          std::cout << "- crossed size: " << csize << std::endl;
          std::cout << "- all_crossed size: " << asize << std::endl;
        }

        const std::size_t num_extra_faces = asize - csize;
        CGAL_assertion(num_extra_faces != 0);

        bool is_ok = true;
        if (num_extra_faces == 2) {

          PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed.back());
          const auto opposite = this->opposite(all_crossed.back(), ivertex);
          const auto& mesh = this->mesh(pvertex);
          auto he = mesh.halfedge(pvertex.second, propagated.second);
          const PEdge qpedge(pvertex.first, mesh.edge(he));
          // std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

          PFace target_pface = null_pface();
          for (const auto pface : pfaces(pvertex.first)) {
            for (const auto pedge : pedges_of_pface(pface)) {
              if (pedge == qpedge) {
                target_pface = pface;
                break;
              }
            }
          }
          CGAL_assertion(target_pface != null_pface());

          const auto tt = pedges_of_pface(target_pface);
          std::vector<PEdge> pedges;
          pedges.reserve(tt.size());
          for (const auto t : tt) pedges.push_back(t);

          PEdge other_pedge = null_pedge();
          for (std::size_t j = 0; j < pedges.size(); ++j) {
            if (pedges[j] == qpedge) {
              const std::size_t jp = (j + 1) % pedges.size();
              const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
              const auto& pedge1 = pedges[jm];
              const auto& pedge2 = pedges[jp];
              const auto iv1 = this->ivertex(this->target(pedge1));
              const auto iv2 = this->ivertex(this->source(pedge2));
              if (iv1 == opposite) {
                CGAL_assertion(iv2 != opposite);
                other_pedge = pedge1;
                break;
              } else if (iv2 == opposite) {
                CGAL_assertion(iv1 != opposite);
                other_pedge = pedge2;
                break;
              } else {
                CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
              }
            }
          }
          CGAL_assertion(other_pedge != null_pedge());
          // std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

          IEdge other_iedge;
          const auto& iedges = support_plane(pvertex).iedges();
          CGAL_assertion(has_iedge(other_pedge));
          const auto query_iedge = this->iedge(other_pedge);
          for (const auto& iedge : iedges) {
            if (iedge == query_iedge) continue;
            if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
              if (line_idx(query_iedge) == line_idx(iedge)) {
                other_iedge = iedge;
              }
            }
          }
          CGAL_assertion(other_iedge != null_iedge());
          // std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

          if (!is_occupied(propagated, other_iedge).first) {
            is_ok = false;
          }
        }

        if (is_ok) {
        if (num_extra_faces < 3) {

          CGAL_assertion_msg(false, "TODO: FRONT, CROSSED < LIMIT, 1 or 2 FACES!");
          CGAL_assertion(future_points.size() == asize);
          CGAL_assertion(future_directions.size() == asize);

          for (std::size_t i = csize; i < asize; ++i) {
            if (m_verbose) std::cout << "- adding extra face" << std::endl;

            PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed[i]);
            if (propagated == null_pvertex()) {

              const Line_2 iedge_line = segment_2(pvertex.first, all_crossed[i]).supporting_line();
              const Point_2 pinit = iedge_line.projection(point_2(pvertex));

              const IVertex opposite = this->opposite(all_crossed[i], ivertex);
              Point_2 future_point = to_2d(pvertex.first, opposite);

              Vector_2 future_direction = Vector_2(pinit, future_point);
              future_point = pinit - m_current_time * future_direction;

              // auto tmp = future_direction;
              // tmp = KSR::normalize(tmp);
              // std::cout << "future tmp: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

              const FT dot_product = future_direction * future_directions[i];
              if (dot_product < FT(0)) {
                future_direction = -future_directions[i];
                future_point = pinit - m_current_time * future_direction;
              } else {
                future_direction = future_directions[i];
                future_point = future_points[i];
              }

              // future_point = future_points[i]; // old, does not work
              // future_direction = future_directions[i]; // old, does not work

              propagated = add_pvertex(pvertex.first, future_point);
              direction(propagated) = future_direction;
              new_vertices.push_back(propagated);

              // std::cout << "propagated null: " << point_3(propagated) << std::endl;
              const PFace pface = pface_of_pvertex(pvertex);
              const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
              this->k(new_pface) = this->k(pface);
              previous = propagated;

              const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
              connect(pedge, all_crossed[i]);
              connect(propagated, all_crossed[i]);
              crossed.push_back(all_crossed[i]); // remove events from this one

              CGAL_assertion_msg(false, "TODO: FRONT, NULL PROPAGATED CASE!");

            } else {

              // std::cout << "propagated std: " << point_3(propagated) << std::endl;
              CGAL_assertion(i == asize - 1);

              // Old code!
              // const PFace pface = pface_of_pvertex(pvertex);
              // const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
              // this->k(new_pface) = this->k(pface);
              // previous = propagated;
              // continue;

              // New code!
              const auto opposite = this->opposite(all_crossed[i], ivertex);
              const auto& mesh = this->mesh(pvertex);
              auto he = mesh.halfedge(pvertex.second, propagated.second);
              const PEdge qpedge(pvertex.first, mesh.edge(he));
              // std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

              PFace target_pface = null_pface();
              for (const auto pface : pfaces(pvertex.first)) {
                for (const auto pedge : pedges_of_pface(pface)) {
                  if (pedge == qpedge) {
                    target_pface = pface;
                    break;
                  }
                }
              }
              CGAL_assertion(target_pface != null_pface());

              const auto tt = pedges_of_pface(target_pface);
              std::vector<PEdge> pedges;
              pedges.reserve(tt.size());
              for (const auto t : tt) pedges.push_back(t);

              PEdge other_pedge = null_pedge();
              for (std::size_t j = 0; j < pedges.size(); ++j) {
                if (pedges[j] == qpedge) {
                  const std::size_t jp = (j + 1) % pedges.size();
                  const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
                  const auto& pedge1 = pedges[jm];
                  const auto& pedge2 = pedges[jp];
                  const auto iv1 = this->ivertex(this->target(pedge1));
                  const auto iv2 = this->ivertex(this->source(pedge2));
                  if (iv1 == opposite) {
                    CGAL_assertion(iv2 != opposite);
                    other_pedge = pedge1;
                    break;
                  } else if (iv2 == opposite) {
                    CGAL_assertion(iv1 != opposite);
                    other_pedge = pedge2;
                    break;
                  } else {
                    CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
                  }
                }
              }
              CGAL_assertion(other_pedge != null_pedge());
              // std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

              IEdge other_iedge;
              const auto& iedges = support_plane(pvertex).iedges();
              CGAL_assertion(has_iedge(other_pedge));
              const auto query_iedge = this->iedge(other_pedge);
              for (const auto& iedge : iedges) {
                if (iedge == query_iedge) continue;
                if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
                  if (line_idx(query_iedge) == line_idx(iedge)) {
                    other_iedge = iedge;
                  }
                }
              }
              CGAL_assertion(other_iedge != null_iedge());
              // std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

              // if (!is_occupied(propagated, other_iedge).first) {
              //   break;
              // }

              CGAL_assertion(m_points.find(std::make_pair(pvertex.first, other_iedge)) != m_points.end());
              CGAL_assertion(m_directions.find(std::make_pair(pvertex.first, other_iedge)) != m_directions.end());
              Point_2 future_point = m_points.at(std::make_pair(pvertex.first, other_iedge));
              Vector_2 future_direction = m_directions.at(std::make_pair(pvertex.first, other_iedge));

              // const Point_2 pinit = point_2(propagated);
              // Point_2 future_point = point_2(pvertex.first, this->opposite(other_iedge, opposite));
              // Vector_2 future_direction = Vector_2(pinit, future_point);
              // future_point = pinit - m_current_time * future_direction;

              auto tmp = future_direction;
              tmp = KSR::normalize(tmp);
              // std::cout << "future tmp: " << to_3d(pvertex.first, point_2(propagated) + m_current_time * tmp) << std::endl;

              const auto before = propagated;
              propagated = add_pvertex(propagated.first, future_point);
              direction(propagated) = future_direction;


              // std::cout << "before: " << point_3(before) << std::endl;
              // std::cout << "propagated: " << point_3(propagated) << std::endl;

              std::size_t count = 0;
              for (const IVertex& iver : { this->source(other_iedge), this->target(other_iedge) }) {
                const Point_2 pi = to_2d(propagated.first, iver);
                const Segment_2 sv(
                  point_2(propagated, m_current_time),
                  point_2(propagated, max_time));
                if (sv.to_vector() * Vector_2(sv.source(), pi) < FT(0)) {
                  ++count;
                  continue;
                }
              }
              if (count == 2) {

                const PFace pface = pface_of_pvertex(pvertex);
                const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, before});
                this->k(new_pface) = this->k(pface);
                previous = before;
                support_plane(pvertex.first).remove_vertex(propagated.second);
                // break;

                // const Point_2 pinit = point_2(before);
                // future_point = point_2(pvertex.first, this->opposite(other_iedge, opposite));
                // future_direction = Vector_2(pinit, future_point);
                // future_point = pinit - m_current_time * future_direction;

                // support_plane(propagated).set_point(propagated.second, future_point);
                // direction(propagated) = future_direction;

                CGAL_assertion_msg(false, "TODO! DOES IT WORK AT ALL?");
              } else {
                new_vertices.push_back(propagated);
              }

              const PFace pface = pface_of_pvertex(pvertex);
              const PFace new_pface = add_pface(std::array<PVertex, 4>{pvertex, previous, propagated, before});
              this->k(new_pface) = this->k(pface);
              previous = propagated;

              CGAL_assertion_msg(false, "DEBUG THIS CASE!");

              const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(before.second, propagated.second));
              connect(pedge, other_iedge);
              connect(propagated, other_iedge);
              crossed.push_back(other_iedge);
            }
          }
          CGAL_assertion_msg(false, "TODO: TEST THIS LOOP!");
        } else {

          // std::cout << "crossed size: " << crossed.size() << std::endl;
          // std::cout << "all crossed size: " << all_crossed.size() << std::endl;
          // for (const auto& iedge : all_crossed) {
          //   std::cout << segment_3(iedge) << std::endl;
          // }

          CGAL_assertion_msg(false, "TODO: FRONT, CROSSED < LIMIT, MULTIPLE FACES!");
        }
        } else {

        }
      }

      if (crossed.size() == all_crossed.size()) {
        // continue...
      }

      if (crossed.size() > all_crossed.size()) {
        // continue ..
        // std::cout << "crossed size: " << crossed.size() << std::endl;
        // std::cout << "all crossed size: " << all_crossed.size() << std::endl;
        CGAL_assertion_msg(false, "TODO: FRONT, CROSSED > LIMIT!");
      }
    }
    else // Open case
    {
      if (m_verbose) {
        std::cout << "*** OPEN CASE" << std::endl;
      }

      // const Direction_2 dir_prev(point_2(prev) - point_2(pvertex));
      // const Direction_2 dir_next(point_2(next) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      const FT next_time = last_event_time(next);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));
      CGAL_assertion(next_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 dir_prev(tmp_prev - point_2(pvertex.first, ivertex));
      const Direction_2 dir_next(tmp_next - point_2(pvertex.first, ivertex));

      // std::cout << to_3d(prev.first, tmp_prev) << std::endl;
      // std::cout << to_3d(next.first, tmp_next) << std::endl;

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (dir_next.counterclockwise_in_between(
          iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;
        const Direction_2& dir = iedges[iedge_idx].second;

        if (!dir.counterclockwise_in_between (dir_next, dir_prev))
          break;

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);

        iedge_idx = (iedge_idx + 1) % iedges.size();

        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: OPEN WHY SO MANY ITERATIONS?");
        }
        ++iter;
      }

      CGAL_assertion(crossed.size() != 0);

      if (m_verbose) {
        std::cout << "- crossed " << crossed.size() << " iedges: " << std::endl;
        for (const auto& iedge : crossed) {
          std::cout << segment_3(iedge) << std::endl;
        }
      }

      std::vector<Point_2> future_points(crossed.size());
      std::vector<Vector_2> future_directions(crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        const bool is_parallel = compute_future_point_and_direction(
          pvertex, prev, next, crossed[i], future_points[i], future_directions[i]);

        if (is_parallel) {
          if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
            prev_iedge = crossed[i];
          }
          if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
            next_iedge = crossed[i];
          }
        }
      }
      CGAL_assertion(future_points.size() == crossed.size());
      CGAL_assertion(future_directions.size() == crossed.size());

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "open saved: " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            pvertex, prev, next, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      {
        PVertex cropped;
        if (next_iedge != null_iedge() && next_iedge == crossed.front()) {
          if (m_verbose) std::cout << "- next parallel case" << std::endl;

          cropped = next;
          Point_2 future_point; Vector_2 future_direction;
          const auto pair = this->border_prev_and_next(next);
          const auto nnext = pair.second;
          compute_future_point_and_direction(
            0, next, nnext, next_iedge, future_point, future_direction);
          future_points[0] = future_point;
          future_directions[0] = future_direction;

        } else {
          if (m_verbose) std::cout << "- standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
          // future_point = future_points.front();
          // future_direction = future_directions.front();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.front());
        connect(cropped, crossed.front());

        support_plane(cropped).set_point(cropped.second, future_points.front());
        direction(cropped) = future_directions.front();

        if (m_verbose) {
          // std::cout << "direction cropped 1: " << direction(cropped) << std::endl;
          std::cout << "- cropped 1: " << point_3(cropped) << std::endl;
        }
      }

      for (std::size_t i = 1; i < crossed.size() - 1; ++i)
      {
        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        connect(propagated, crossed[i]);
        new_vertices.push_back(propagated);
        if (m_verbose) {
          std::cout << "- propagated " << std::to_string(i) << ": " << point_3(propagated) << std::endl;
        }
      }

      {
        PVertex cropped;
        if (prev_iedge != null_iedge() && prev_iedge == crossed.back()) {
          if (m_verbose) std::cout << "- prev parallel case" << std::endl;

          cropped = prev;
          Point_2 future_point; Vector_2 future_direction;
          const auto pair = this->border_prev_and_next(prev);
          const auto pprev = pair.first;
          compute_future_point_and_direction(
            0, prev, pprev, prev_iedge, future_point, future_direction);
          future_points[future_points.size() - 1] = future_point;
          future_directions[future_directions.size() - 1] = future_direction;

        } else {
          if (m_verbose) std::cout << "- standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
          // future_point = future_points.back();
          // future_direction = future_directions.back();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.back());
        connect(cropped, crossed.back());

        support_plane(cropped).set_point(cropped.second, future_points.back());
        direction(cropped) = future_directions.back();

        if (m_verbose) {
          // std::cout << "direction cropped 2: " << direction(cropped) << std::endl;
          std::cout << "- cropped 2: " << point_3(cropped) << std::endl;
        }
      }

      if (m_verbose) std::cout << "- new pvertices size: " << new_vertices.size() << std::endl;
      CGAL_assertion(new_vertices.size() == crossed.size());

      bool is_occupied_edge_back, bbox_reached_back;
      std::tie(is_occupied_edge_back, bbox_reached_back) = is_occupied(pvertex, ivertex, crossed.back());
      // std::tie(is_occupied_edge_back, bbox_reached_back) = collision_occured(pvertex, crossed.back());

      if (m_verbose) {
        std::cout << "- is already occupied back / bbox: " << is_occupied_edge_back << "/" << bbox_reached_back << std::endl;
      }

      bool is_occupied_edge_front, bbox_reached_front;
      std::tie(is_occupied_edge_front, bbox_reached_front) = is_occupied(pvertex, ivertex, crossed.front());
      // std::tie(is_occupied_edge_front, bbox_reached_front) = collision_occured(pvertex, crossed.front());

      if (m_verbose) {
        std::cout << "- is already occupied fron / bbox: " << is_occupied_edge_front << "/" << bbox_reached_front << std::endl;
      }

      const auto pface = pface_of_pvertex(pvertex);
      if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
      if (bbox_reached_back) {

        CGAL_assertion(bbox_reached_front);
        if (m_verbose) std::cout << "- stop bbox back" << std::endl;

      } else if (bbox_reached_front) {

        CGAL_assertion(bbox_reached_back);
        if (m_verbose) std::cout << "- stop bbox front" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) == 1) {

        if (m_verbose) std::cout << "- stop back && front k = 1" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) > 1) {

        this->k(pface)--;
        CGAL_assertion(this->k(pface) >= 1);
        add_new_pfaces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        if (m_verbose) std::cout << "- continue back && front k > 1" << std::endl;

      } else if ((!is_occupied_edge_back && !is_occupied_edge_front)) {

        add_new_pfaces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        if (m_verbose) std::cout << "- continue !back && !front" << std::endl;

      } else if (is_occupied_edge_back || is_occupied_edge_front) {

        // if (this->k(pface) > 1) {
        //   this->k(pface)--;
        // }
        // CGAL_assertion(this->k(pface) >= 1);
        add_new_pfaces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        if (m_verbose) std::cout << "- continue back || front" << std::endl;

        // std::cout << "pv pface: "   << str(pface_of_pvertex(pvertex))      << std::endl;
        // std::cout << "back pface: " << str(pface_of_pvertex(pvertices[1])) << std::endl;
        // std::cout << "fron pface: " << str(pface_of_pvertex(pvertices[2])) << std::endl;
        // CGAL_assertion_msg(false, "TEST THIS CASE: BACK || FRONT!");

      } else {
        CGAL_assertion_msg(false, "TODO: ADD NEW OPEN CASE! DO NOT FORGET TO UPDATE K!");
      }

      if (m_verbose) {
        // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
        std::cout << "- k intersections after: " << this->k(pface) << std::endl;
      }

      // for (std::size_t i = 1; i < crossed.size() - 1; ++i) {
      //   PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, new_vertices[i].second));
      //   connect(pedge, crossed[i]);
      //   connect(new_vertices[i], crossed[i]);
      // }
    }

    support_plane(support_plane_idx).remove_vertex(front.second);
    support_plane(support_plane_idx).remove_vertex(back.second);

    // push also remaining vertex so that its events are recomputed
    // std::cout << "pushing new pv: " << str(pvertex) << std::endl;
    // std::cout << "pv direction: " << direction(pvertex) << std::endl;
    new_vertices.push_back(pvertex);
    crossed.push_back(iedge(pvertex));

    if (m_verbose) {
      std::cout << "- new pvertices:";
      for (const PVertex& pv : new_vertices)
        std::cout << " " << str(pv);
      std::cout << std::endl;
    }

    // if (has_iedge(prev) && !is_frozen(prev)) {
    // // if (iedge(prev) != iedge(pvertex)) {
    //   std::cout << "pushing new prev: " << str(prev) << std::endl;
    //   new_vertices.push_back (prev);
    // }

    // if (has_iedge(next) && !is_frozen(next)) {
    // // if (back_constrained) {
    //   std::cout << "pushing new next: " << str(next) << std::endl;
    //   new_vertices.push_back (next);
    // }

    return new_vertices;
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_EXPERIMENTAL_H
