// Copyright (c) 2009  INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Pierre Alliez, Camille Wormser

#ifndef CGAL_AABB_SEARCH_TREE_H
#define CGAL_AABB_SEARCH_TREE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>

namespace CGAL
{

template <class Traits>
struct AABB_search_tree
{

  typedef typename Traits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename Point_and_primitive_id::first_type Point;
  typedef typename Point_and_primitive_id::second_type Id;
  typedef First_of_pair_property_map<Point_and_primitive_id> Pmap;
  typedef Search_traits_adapter<Point_and_primitive_id, Pmap, Traits> TreeTraits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
private:
  Tree m_tree;

  Point_and_primitive_id get_p_and_p(const Point_and_primitive_id& p)
  {
    return p;
  }

  Point_and_primitive_id get_p_and_p(const Point& p)
  {
    return Point_and_primitive_id(p, Id());
  }

public:
  template <class ConstPointIterator>
  AABB_search_tree(ConstPointIterator begin, ConstPointIterator beyond)
      : m_tree{}
  {
    std::vector<Point_and_primitive_id> points;
    while (begin != beyond) {
      Point_and_primitive_id pp = get_p_and_p(*begin);
      points.emplace_back(pp);
      ++begin;
    }
    m_tree.insert(points.begin(), points.end());
    m_tree.build();
  }

  template <typename Point>
  Point_and_primitive_id closest_point(const Point& query) const
  {
    Neighbor_search search(m_tree, query, 1);
    return search.begin()->first;
  }
};

}

#endif // CGAL_AABB_SEARCH_TREE_H
