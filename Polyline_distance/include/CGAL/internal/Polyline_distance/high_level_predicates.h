#pragma once

#include "low_level_predicates.h"
#include "predicate_types.h"

#include <iterator>

namespace HLPred
{

Interval intersection_interval(const Point& circle_center, distance_t radius, Point line_start, Point line_end)
{
	using OutputType = std::pair<Circular_arc_point_2, unsigned>;

	auto circle_2 = Circle(circle_center, radius);
	auto line_arc_2 = LineArc(line_start, line_end);

	std::vector<OutputType> intersections;
	LLPred::IntersectionAlgorithm::intersection(circle_2, line_arc_2, std::back_inserter(intersections));

	std::vector<distance_t> ratios;
	for (auto const& intersection: intersections) {
		assert(line_start.x != line_end.x || line_start.y != line_end.y);

		distance_t ratio;
		if (line_start.x != line_end.x) {
			ratio = (intersection.first.x - line_start.x)/(line_end.x - line_start.x);
		}
		else {
			ratio = (intersection.first.y - line_start.y)/(line_end.y - line_start.y);
		}
		ratios.push_back(ratio);
	}

	assert(intersections.size() <= 2);

	switch (intersections.size())
	{
	case 0:
		return Interval(); // empty interval
	case 1:
		if (intersections[0].second == 2) { return Interval(ratios[0], ratios[0]); }
		// TODO: replace by predicates
		if (line_start.dist_sqr(circle_center) <= radius*radius) { return Interval(0., ratios[0]); }
		if (line_end.dist_sqr(circle_center) <= radius*radius) { return Interval(ratios[0], 1.); }
		assert(false);
	case 2:
		return Interval(std::min(ratios[0], ratios[1]), std::max(ratios[0], ratios[1]));
	}
}

} // end HLPred namespace
