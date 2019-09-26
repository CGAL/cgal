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

#ifndef CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
#define CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H

#include <CGAL/basic.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Dimension.h>

// FIXME: vector too restricted?
#include <iterator>
#include <vector>

namespace CGAL{

template <typename PointRange>
using PointRangeKernel = typename CGAL::Kernel_traits<
                           typename std::iterator_traits<
                             typename PointRange::iterator>::value_type>::Kernel;

template <class PointRange,
          class Traits = PointRangeKernel<PointRange>>
class FrechetDistanceNearNeighborsDS
{
	using NT = typename Traits::FT;
	using Point = typename Traits::Point_2;
	using Curve = std::vector<Point>;
	using Curves = std::vector<Curve>;
	using CurveIDs = std::vector<std::size_t>; // FIXME: 

public:
	FrechetDistanceNearNeighborsDS() = default;

	void fill(const Curves& curves); // FIXME: should this be a range?
	CurveIDs get_close_curves(const Curve& curve, NT distance);
private:
	Curves curves;

	using TreeTraits = Search_traits_d<Cartesian_d<typename Traits::FT>, Dimension_tag<8>>;
	using TreePoint = typename TreeTraits::Point_d;
	using TreePoints = std::vector<TreePoint>;

	Kd_tree<TreeTraits> kd_tree; // FIXME: what traits?
	TreePoint to_kd_tree_point(const Curve& curve) const;

	struct QueryItem {
		using D = Dimension_tag<8>;
		using Point_d = TreePoint;
		using FT = NT;

		bool contains(Point_d p) const { return true; } // FIXME
		bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
		bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const { return true; } // FIXME
	};
	QueryItem construct_fuzzy_query_item(const Curve& curve, NT distance);
};

// TODO: store preprocessed curves after CGALization
template <class PointRange, class Traits>
void
FrechetDistanceNearNeighborsDS<PointRange,Traits>::fill(const Curves& curves)
{
	for (auto const& curve: curves) {
		kd_tree.insert(to_kd_tree_point(curve));
	}
	kd_tree.build();
}

template <class PointRange, class Traits>
auto
FrechetDistanceNearNeighborsDS<PointRange,Traits>::get_close_curves(const Curve& curve, NT distance)
-> CurveIDs
{
	TreePoints candidates;

	kd_tree.search(std::back_insert_iterator<TreePoints>(candidates), construct_fuzzy_query_item(curve, distance));

	// TODO: convert candidates to curve ids. This is unnecessary but the kd tree doesn't seem to provide the functionality that we want of returning values associated with d-dimensional keys.
	CurveIDs result;

	auto predicate = [&](CurveID id) {
		return !continuous_Frechet_distance_less_than(curve, curves[id], distance);
	};
	auto new_end = std::remove_if(result.begin(), result.end(), predicate);
	result.erase(new_end, result.end());

	return result;
}

template <class PointRange, class Traits>
auto
FrechetDistanceNearNeighborsDS<PointRange,Traits>::to_kd_tree_point(const Curve& curve) const
-> TreePoint
{
	return TreePoint(8, ORIGIN); // FIXME: why is there no point type with dimension as template parameter?
}

template <class PointRange, class Traits>
auto
FrechetDistanceNearNeighborsDS<PointRange,Traits>::construct_fuzzy_query_item(const Curve& curve, NT distance)
-> QueryItem
{
	return QueryItem();
}

} // end of namespace CGAL

#endif // CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
