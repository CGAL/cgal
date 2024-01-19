#pragma once

#include "geometry_basics.h"
#include "predicate_types.h"

#include <CGAL/number_utils.h>
#include <iterator>

distance_t distance(Point const& p, Point const& q)
{
	return CGAL::sqrt(CGAL::squared_distance(p, q));
}

namespace LLPred
{

//
// compare_squared_distance
//

enum Comparison_result {
	LESS,
	EQUAL,
	LARGER,
};

Comparison_result compare_squared_distance(Point const& p, Point const& q, distance_t d2)
{
	auto const pq_dist_sqr = CGAL::squared_distance(p, q);

	if (pq_dist_sqr < d2) { return LESS; }
	if (pq_dist_sqr == d2) { return EQUAL; }
	else { return LARGER; }
}

//
// circle-segment intersection
//

using IntersectionOutputType = std::pair<Circular_arc_point_2, unsigned>;
using IntersectionBackInserter = std::back_insert_iterator<std::vector<IntersectionOutputType>>;

class IntersectionAlgorithm
{
	static constexpr distance_t eps = 1e-8;

public:
   /*
    * Returns which section of the line segment from line_start to line_end is inside the circle given by circle_center and radius.
    * If the circle and line segment do not intersect, the result is the empty Interval (and outer is the empty Interval, too).
	* Otherwise the result is an interval [x,y] such that the distance at x and at y is at most the radius, i.e., [x,y] is a subset of the free interval. 
	* The optional output "outer" is an interval strictly containing the free interval.
	* In other words, "outer" is an interval [x',y'] containing [x,y] such that x-x', y'-y <= eps and:
	* If x = 0 then x' = -eps, while if x > 0 then the distance at x' is more than the radius.
	* If y = 1 then y' = 1+eps, while if y < 1 then the distance at y' is more than the radius.
    */
	static void intersection(Circle const& circle, LineArc line_arc, IntersectionBackInserter it, Interval* outer);
private:
	IntersectionAlgorithm() {} // Make class static-only
	static inline bool smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr);
	static inline distance_t distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center);

	static constexpr distance_t save_eps = 0.5 * eps;
	static constexpr distance_t save_eps_half = 0.25 * eps;
};

inline bool IntersectionAlgorithm::smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr) {
	return toOldPoint(circle_center).dist_sqr(toOldPoint(line_start) * (1. - interpolate) + toOldPoint(line_end) * interpolate) <= radius_sqr;
}

inline distance_t IntersectionAlgorithm::distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center) {
	return toOldPoint(circle_center).dist_sqr(toOldPoint(line_start) * (1. - interpolate) + toOldPoint(line_end) * interpolate);
}


void IntersectionAlgorithm::intersection(Circle const& circle, LineArc line_arc, IntersectionBackInserter it, Interval* outer)
{
	Point const& circle_center = circle.center;
	distance_t radius = circle.radius;
	Point const& line_start = line_arc.start;
	Point const& line_end = line_arc.end;

    // The line can be represented as line_start + lambda * v
    const OldPoint v = toOldPoint(line_end) - toOldPoint(line_start);
	const distance_t rad_sqr = radius * radius;
	
    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x()^2 + p.y()^2) = radius
    // <=> p.x()^2 + p.y()^2 = radius^2
    // <=> (line_start.x() + lambda * v.x())^2 + (line_start.y() + lambda * v.y())^2 = radius^2
    // <=> (line_start.x()^2 + 2 * line_start.x() * lambda * v.x() + lambda^2 * v.x()^2) + (line_start.y()^2 + 2 * line_start.y() * lambda * v.y() + lambda^2 * v.y()^2) = radius^2
    // <=> lambda^2 * (v.x()^2 + v.y()^2) + lambda * (2 * line_start.x() * v.x() + 2 * line_start.y() * v.y()) + line_start.x()^2 + line_start.y()^2) - radius^2 = 0
    // let a := v.x()^2 + v.y()^2, 
	// let b := line_start.x() * v.x() + line_start.y() * v.y(), 
	// let c := line_start.x()^2 + line_start.y()^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
	
    const distance_t a = pow2(v.x()) + pow2(v.y());
    const distance_t b = (line_start.x() - circle_center.x()) * v.x() + (line_start.y() - circle_center.y()) * v.y();
    const distance_t c = pow2(line_start.x() - circle_center.x()) + pow2(line_start.y() - circle_center.y()) - pow2(radius);

	distance_t mid = - b / a;
    distance_t discriminant = pow2(mid) - c / a;

	const bool smallDistAtZero = smallDistanceAt(0., line_start, line_end, circle_center, rad_sqr);
	const bool smallDistAtOne = smallDistanceAt(1., line_start, line_end, circle_center, rad_sqr);
	bool smallDistAtMid = smallDistanceAt(mid, line_start, line_end, circle_center, rad_sqr);
	
	if (smallDistAtZero && smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(-eps, 1. + eps); }

		it = std::make_pair(line_start, 1);
		it = std::make_pair(line_end, 1);
		return;
	}
	
	if (!smallDistAtMid && smallDistAtZero) {
		mid = 0.;
		smallDistAtMid = true;
	}
	else if (!smallDistAtMid && smallDistAtOne) {
		mid = 1.;
		smallDistAtMid = true;
	}
	
	// Here we need the guarantee that if the free interval has length at least eps
	// then at mid the distance is <=radius
	// This is an assumption about the precision of distance_t computations
	// All remaining rules are free of such assumptions! 
	// (except for trivial ones like this: x + y and x - y have distance at most 2y up to negligible error)
    if (!smallDistAtMid) {
		if (outer != nullptr) { *outer = Interval(); }
		return; // no intersection;
    }
	
	if (mid <= 0. and !smallDistAtZero) {
		if (outer != nullptr) { *outer = Interval(); }
		return;
	}
	if (mid >= 1. and !smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(); }
		return;
	}
	
	discriminant = std::max<distance_t>(discriminant, 0.);
	distance_t sqrt_discr = 0.;
	bool sqrt_discr_computed = false;
	distance_t begin, end;
	
	if (smallDistAtZero) {
		begin = 0.;
		if (outer != nullptr) { outer->begin = -eps; }
	}
	else {
		sqrt_discr = std::sqrt(discriminant);
		sqrt_discr_computed = true;
		
		const distance_t lambda1 = mid - sqrt_discr;
		const distance_t innershift = std::min<distance_t>(lambda1 + save_eps_half, std::min<distance_t>(1., mid));
		const distance_t outershift = lambda1 - save_eps_half;
		if (innershift >= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			begin = innershift;
			if (outer != nullptr) { outer->begin = outershift; }
		}
		else {
			distance_t left = 0., right = std::min<distance_t>(mid, 1.);
			// invariants throughout binary search:
			//  * !smallDistanceAt(left)
			//  * smallDistanceAt(right)
			//  * 0 <= left <= right <= min(mid,1)
			// Clearly this is stays true after an iteration.
			// Why is it true in the beginning?
			// If smallDistanceAt(0.) then begin would already be set (fourth rule).
			// If !smallDistanceAt(right), then either !smallDistanceAt(mid), contradicting the very first rule, 
			//  or mid >= 1. and smallDistanceAt(1.), contradicting the third rule.
			// Finally, since !smallDistanceAt(left) we cannot have mid <= 0 by the second rule. Thus, right = min(mid,1) >= 0. = left
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { right = m; }
				else { left = m; }
			}
			begin = right;
			if (outer != nullptr) { outer->begin = left; }
		}
	}
	
	if (smallDistAtOne) {
		end = 1.;
		if (outer != nullptr) { outer->end = 1. + eps; }
	}
	else {
		if (!sqrt_discr_computed) {
			sqrt_discr = std::sqrt(discriminant);
		}
		
		const distance_t lambda2 = mid + sqrt_discr;
		const distance_t innershift = std::max<distance_t>(lambda2 - save_eps_half, std::max<distance_t>(0., mid));
		const distance_t outershift = lambda2 + save_eps_half;
		if (innershift <= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			end = innershift;
			if (outer != nullptr) { outer->end = outershift; }
		}
		else {
			distance_t left = std::max<distance_t>(mid, 0.), right = 1.;
			// invariants throughout binary search:
			//  * smallDistanceAt(left)
			//  * !smallDistanceAt(right)
			//  * max(mid,0) <= left <= right <= 1
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { left = m; }
				else { right = m; }
			}
			end = left;
			if (outer != nullptr) { outer->end = right; }
		}
	}
	
	assert(smallDistanceAt(begin, line_start, line_end, circle_center, rad_sqr));
	assert(smallDistanceAt(end, line_start, line_end, circle_center, rad_sqr));
	assert(0. <= begin && begin <= end && end <= 1.);
	
	assert(outer == nullptr || outer->begin < 0. || !smallDistanceAt(outer->begin, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || outer->end > 1. || !smallDistanceAt(outer->end, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || (outer->begin < begin && begin - outer->begin <= eps));
	assert(outer == nullptr || (outer->end > end && outer->end - end <= eps));
	
	assert(outer == nullptr || (outer->begin <= 1.));
	assert(outer == nullptr || (outer->end >= 0.));
	
	// These asssertions fail - use exact arithmetic to make it work??
	//assert(begin - eps < 0. || !smallDistanceAt(begin - eps, line_start, line_end, circle_center, rad_sqr));
	//assert(end + eps > 1. || !smallDistanceAt(end + eps, line_start, line_end, circle_center, rad_sqr));
	
	if (begin == end) {
		OldPoint old_p = toOldPoint(line_start) + (toOldPoint(line_end) - toOldPoint(line_start))*begin;
		it = std::make_pair(Point(old_p.x(), old_p.y()), 2);
	}
	else {
		OldPoint old_p1 = toOldPoint(line_start) + (toOldPoint(line_end) - toOldPoint(line_start))*begin;
		it = std::make_pair(Point(old_p1.x(), old_p1.y()), 1);
		OldPoint old_p2 = toOldPoint(line_start) + (toOldPoint(line_end) - toOldPoint(line_start))*end;
		it = std::make_pair(Point(old_p2.x(), old_p2.y()), 1);
	}
}

void intersection(Circle const& circle, LineArc line_arc, IntersectionBackInserter it, Interval* outer)
{
	IntersectionAlgorithm::intersection(circle, line_arc, it, outer);
}

} // end LLPred namespace
