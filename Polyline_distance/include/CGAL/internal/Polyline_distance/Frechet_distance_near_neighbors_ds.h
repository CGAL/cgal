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
	Point_d to_kd_tree_point(const Curve& curve) const;

	class QueryItem {
		const Curve& curve;
		const NT distance;

	public:
		using D = FrechetKdTree::D;
		using Point_d = FrechetKdTree::Point_and_id;
		using FT = NT;

		QueryItem(Curve const& curve, NT distance) : curve(curve), distance(distance) {}

		bool contains(Point_d p) const { return true; } // FIXME
		bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
		bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
	};
};

template <class Traits>
auto
FrechetKdTree<Traits>::to_kd_tree_point(const Curve& curve) const
-> Point_d
{
	CGAL_precondition(!curve.empty());

	// TODO: remove this when curves are preprocessed first
	NT x_min, y_min, x_max, y_max;
	x_min = y_min = std::numeric_limits<NT>::max();
	x_max = y_max = std::numeric_limits<NT>::min();
	for (auto const& point: curve) {
		x_min = std::min(x_min, point.x());
		y_min = std::min(y_min, point.y());
		x_max = std::max(x_max, point.x());
		y_max = std::max(y_max, point.y());
	}

	// FIXME: this is ugly......
	std::vector<NT> values = {
		curve.front().x(),
		curve.front().y(),
		curve.back().x(),
		curve.back().y(),
		x_min,
		y_min,
		x_max,
		y_max
	};
	return Point_d(D::value, values.begin(), values.end());
}

template <class Traits>
void
FrechetKdTree<Traits>::insert(Curves const& curves)
{
	for (std::size_t id = 0; id < curves.size(); ++id) {
		auto kd_tree_point = to_kd_tree_point(curves[id]);
		kd_tree.insert(Point_and_id{kd_tree_point, id});
	}
}

template <class Traits>
void
FrechetKdTree<Traits>::build()
{
	kd_tree.build();
}

template <class Traits>
CurveIDs
FrechetKdTree<Traits>::search(Curve const& curve, NT distance)
{
	using Candidates = std::vector<Point_and_id>;
	Candidates candidates;

	kd_tree.search(std::back_insert_iterator<Candidates>(candidates), QueryItem(curve, distance));

	// FIXME: ugly.........
	CurveIDs result;
	for (auto const& candidate: candidates) {
		result.push_back(get<1>(candidate));
	}
	return result;
}

} // end of namespace CGAL

#endif // CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
