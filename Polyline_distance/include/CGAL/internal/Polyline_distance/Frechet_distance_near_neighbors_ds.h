// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
#define CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H

#include <CGAL/basic.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Dimension.h>

// FIXME: vector too restricted?
#include <iterator>
#include <vector>

#include <boost/iterator/zip_iterator.hpp>

namespace CGAL
{

// TODO: which boost parts can be replaced by modern C++ std library?
template <class Traits>
class FrechetKdTree
{
	// FIXME: somehow make it more clear that those come from the Fréchet things...
	using NT = typename Traits::FT;
	using Point = typename Traits::Point_2;
	using Curve = std::vector<Point>;
	using Curves = std::vector<Curve>;
	using CurveIDs = std::vector<std::size_t>; // FIXME: 

	using D = Dimension_tag<8>;
	using Tree_traits_base = Search_traits_d<Cartesian_d<typename Traits::FT>, D>;
	using Point_d = typename Tree_traits_base::Point_d;
	using Point_and_id = boost::tuple<Point_d,std::size_t>; // TODO: use some ID type here instead of size_t
	using Tree_traits = CGAL::Search_traits_adapter<Point_and_id, CGAL::Nth_of_tuple_property_map<0, Point_and_id>, Tree_traits_base>;
	using Tree = Kd_tree<Tree_traits>;

public:
	FrechetKdTree() = default;

	void insert(Curves const& curves); // TODO: also add other methods of inserting points
	void build();
	CurveIDs search(Curve const& curve, NT distance);

private:
	Kd_tree<Tree_traits> kd_tree; // FIXME: what traits?

	// FIXME: rename
	Point_and_id to_kd_tree_point(const Curve& curve) const;

	struct QueryItem {
		using D = FrechetKdTree::D;
		using Point_d = FrechetKdTree::Point_d;
		using FT = NT;

		bool contains(Point_d p) const { return true; } // FIXME
		bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
		bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
	};
	QueryItem construct_fuzzy_query_item(const Curve& curve, NT distance);
};

template <class Traits>
auto
FrechetKdTree<Traits>::to_kd_tree_point(const Curve& curve) const
-> Point_and_id
{
	return Point_and_id{Point_d(), 0}; // FIXME
}

template <class Traits>
void
FrechetKdTree<Traits>::insert(Curves const& curves)
{
	for (std::size_t id = 0; id < curves.size(); ++id) {
		auto const& curve = curves[id];
		auto kd_tree_point = to_kd_tree_point(curve);
		kd_tree.insert(kd_tree_point);
	}
}

template <class Traits>
void
FrechetKdTree<Traits>::build()
{
	// TODO
}

template <class Traits>
CurveIDs
FrechetKdTree<Traits>::search(Curve const& curve, NT distance)
{
	// TODO
	// kd_tree.search(std::back_insert_iterator<TreePoints>(candidates), construct_fuzzy_query_item(curve, distance));
}

template <class Traits>
auto
FrechetKdTree<Traits>::construct_fuzzy_query_item(const Curve& curve, NT distance)
-> QueryItem
{
	// TODO
	return QueryItem();
}


} // end of namespace CGAL

#endif // CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
