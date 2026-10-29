// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_REGROUP_POLYLINES_H
#define CGAL_REGROUP_POLYLINES_H

#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {

template<class FT>
class Box_for_polylines : public Box_intersection_d::Box_d< FT, 2, Box_intersection_d::ID_NONE>
{
    std::size_t m_gid, m_pid, m_sid;
public:
    typedef Box_intersection_d::Box_d< FT, 2, Box_intersection_d::ID_NONE> Base;
    typedef FT                            NT;
    typedef std::size_t                   ID;

    Box_for_polylines(const Bbox_2& b, std::size_t gid, std::size_t pid, std::size_t sid)
        : Base(b)
        , m_gid(gid)
        , m_pid(pid)
        , m_sid(sid)
    {}

    ID  id() const { return m_gid; }
    ID  polyline_id() const { return m_pid; }
    ID  segment_id() const { return m_sid; }
};

template <class GT>
void regroup_poly(const std::vector<std::vector<typename GT::Point_2> >& polylines,
                  std::vector< std::vector<std::size_t> >& groups,
                  std::vector<int>& nesting_level,
                  std::vector<bool>& does_polyline_self_intersect,
                  std::vector<bool>& is_polygon)
{
  using Point_2 = typename GT::Point_2;
  using Segment_2 = typename GT::Segment_2;

  std::size_t nb_poly = polylines.size();

  // create the bbox for all segments
  std::size_t i = 0;
  std::vector< Box_for_polylines<double> > boxes; // TODO add reserve
  Bbox_2 gbox;

  typename GT::Construct_bbox_2 bbox_2;
  for (std::size_t pid = 0; pid < nb_poly; ++pid)
  {
    const std::vector<Point_2>& polyline = polylines[pid];
    std::size_t ns = polyline.size();
    for (std::size_t k = 0; k < ns - 1; ++k)
    {
      Bbox_2 b2 = bbox_2(polyline[k]) + bbox_2(polyline[k + 1]);
      gbox += b2;
      boxes.push_back(Box_for_polylines<double>(b2, i++, pid, k));
    }
  }

  does_polyline_self_intersect.assign(nb_poly, false);
  is_polygon.resize(nb_poly);
  std::set< std::pair<std::size_t, std::size_t> > intersecting_polylines; // TODO (add as output?)

  for (std::size_t pid = 0; pid < nb_poly; ++pid)
  {
    is_polygon[pid] = polylines[pid].front() == polylines[pid].back();
  }

  auto segment_2 = [&](std::size_t pid, std::size_t sid)
  {
    const Point_2& s = polylines[pid][sid];
    const Point_2& t = polylines[pid][sid + 1];
    return Segment_2(s, t);
  };

  auto lambda = [&](const Box_for_polylines<double>& a, const Box_for_polylines<double>& b)
  {
    Segment_2 sa = segment_2(a.polyline_id(), a.segment_id());
    Segment_2 sb = segment_2(b.polyline_id(), b.segment_id());

    auto inter = typename GT::Intersect_2()(sa, sb);

    if (inter != boost::none)
    {
      //TODO: use internal function, maybe even do_intersect is enough
      if (/* const Point_2* pt =  */boost::get<Point_2>(&*inter))
      {
        if (a.polyline_id() == b.polyline_id())
        {
            if ((a.segment_id() + 1) % (polylines[a.polyline_id()].size() - 1) != b.segment_id() &&
                (b.segment_id() + 1) % (polylines[a.polyline_id()].size() - 1) != a.segment_id())
            {
              does_polyline_self_intersect[a.polyline_id()] = true;
            }
        }
        else
        {
            intersecting_polylines.insert(make_sorted_pair(a.polyline_id(), b.polyline_id()));
        }
      }
      else
        if (/* const Segment_2* s =  */boost::get<Segment_2>(&*inter))
        {
          if (a.polyline_id() == b.polyline_id())
            does_polyline_self_intersect[a.polyline_id()] = true;
          else
            intersecting_polylines.insert(make_sorted_pair(a.polyline_id(), b.polyline_id()));
        }
    }
  };

  // Run the self intersection algorithm with all defaults
  box_self_intersection_d(boxes.begin(), boxes.end(), lambda);

#ifdef CGAL_REGROUP_POLYLINES_VERBOSE
  std::cout << "Self-intersecting polylines:";
  for (std::size_t pid = 0; pid < nb_poly; ++pid)
    if (does_polyline_self_intersect[pid])
      std::cout << " " << pid;
  std::cout << "\n";

  std::cout << "Intersection between polylines:\n";
  for (const std::pair<std::size_t, std::size_t>& p : intersecting_polylines)
    std::cout << "   - (" << p.first << "," << p.second << ")\n";

  std::cout << "Polygons:\n";
  for (std::size_t pid = 0; pid < nb_poly; ++pid)
    if (is_polygon[pid])
      std::cout << " " << pid;
  std::cout << "\n";
#endif

  // min point per polygon
  typename GT::Less_xy_2 less;

  Point_2 max_pt(gbox.xmax(), gbox.ymax(), 0);
  std::vector<Point_2> min_point(nb_poly, max_pt);
  for (std::size_t pid = 0; pid < nb_poly; ++pid)
  {
    for (const Point_2& pt : polylines[pid])
      if (less(pt, min_point[pid]))
        min_point[pid] = pt;
  }

  std::vector<bool> handled = does_polyline_self_intersect;
  // skip also polylines for now
  for (std::size_t pid = 0; pid < nb_poly; ++pid)
    if (!is_polygon[pid])
      handled[pid] = true;

  auto get_extreme_id = [&polylines, nb_poly, less, &min_point](const std::vector<bool>& handled)
  {
    std::size_t xtrm_id(-1);
    for (std::size_t pid = 0; pid < nb_poly; ++pid)
    {
      if (!handled[pid])
      {
        if (xtrm_id == std::size_t(-1))
          xtrm_id = pid;
        else
            if (less(min_point[pid], min_point[xtrm_id]))
              xtrm_id = pid;
      }
    }
    return xtrm_id;
  };

  nesting_level.assign(nb_poly, -1);
  while (true)
  {
    std::size_t xtrm_id = get_extreme_id(handled);

    if (xtrm_id == std::size_t(-1)) break;
    groups.resize(groups.size() + 1); // create new group

    handled[xtrm_id] = true;
    groups.back().push_back(xtrm_id);
    nesting_level[xtrm_id] = 0;

    // collect ids inside the current polygon
    std::vector<std::size_t> nested;
    for (std::size_t pid = 0; pid < nb_poly; ++pid)
    {
      if (!handled[pid])
      {
        if (intersecting_polylines.count(make_sorted_pair(xtrm_id, pid)) == 1)
          continue; // considered as disjoint

        if (bounded_side_2(polylines[xtrm_id].begin(),
            std::prev(polylines[xtrm_id].end()),
            polylines[pid].front(),
            GT()) == ON_BOUNDED_SIDE)
        {
          nested.push_back(pid);
          handled[pid] = true;
          groups.back().push_back(pid);
        }
      }
    }

    //compute the nesting_level
    std::vector<std::pair<int, std::vector<std::size_t>> > nesting_stack;
    if (!nested.empty())
    {
      nesting_stack.resize(1);
      nesting_stack.back().first = 1;
      nesting_stack.back().second.swap(nested);
    }

    while (!nesting_stack.empty())
    {
      std::vector<std::size_t> pids;
      int k;
      std::tie(k, pids) = nesting_stack.back();
      nesting_stack.pop_back();

      std::size_t xtrm_id = pids.front();
      for (std::size_t pid : pids)
        if (less(min_point[pid], min_point[xtrm_id]))
          xtrm_id = pid;


      nesting_level[xtrm_id] = k;

      std::vector<std::size_t> k1_nested, k_nested;
      for (std::size_t pid : pids)
      {
        if (pid == xtrm_id) continue;

        if (intersecting_polylines.count(make_sorted_pair(xtrm_id, pid)) == 1)
        {
          k_nested.push_back(pid);
          continue; // considered as disjoint
        }

        if (bounded_side_2(polylines[xtrm_id].begin(),
            std::prev(polylines[xtrm_id].end()),
            polylines[pid].front(),
            GT()) == ON_BOUNDED_SIDE)
        {
          k1_nested.push_back(pid);
        }
        else
          k_nested.push_back(pid);
      }

      if (!k1_nested.empty())
        nesting_stack.emplace_back(k + 1, std::move(k1_nested));
      if (!k_nested.empty())
        nesting_stack.emplace_back(k, std::move(k_nested));
    }
  }

  for (std::size_t pid = 0; pid < nb_poly; ++pid)
  {
    if (!is_polygon[pid] || does_polyline_self_intersect[pid])
    {
      //put in its own group
      groups.emplace_back(1, pid);
      continue;
    }
  }
}

} // end of CGAL namespace

#endif
