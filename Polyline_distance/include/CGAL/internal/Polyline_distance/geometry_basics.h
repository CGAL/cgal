#pragma once

#include "defs.h"
#include "id.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <sstream>

namespace unit_tests { void testGeometricBasics(); }

//
// distance_t
//

using distance_t = double;

//
// Point
//

struct Point {
    distance_t x;
    distance_t y;

	Point() = default;
	Point(distance_t x, distance_t y) : x(x), y(y) {}

	Point& operator-=(const Point& point);
	Point operator-(const Point& point) const;
	Point& operator+=(const Point& point);
	Point operator+(const Point& point) const;
	Point operator*(const distance_t mult) const;
	Point& operator/=(distance_t distance);
	Point operator/(distance_t distance);

	bool operator==(Point const& other) const;
	bool operator!=(Point const& other) const;

	distance_t dist_sqr(const Point& point) const;
	distance_t dist(const Point& point) const;
};
using Points = std::vector<Point>;
using PointID = ID<Point>;

std::ostream& operator<<(std::ostream& out, const Point& p);

struct PointRange {
	PointID begin;
	PointID end;
};

struct ContinuousPoint {
	PointID point;
	distance_t fraction;

	bool operator<(ContinuousPoint const& other) const {
		return point < other.point || (point == other.point && fraction < other.fraction);
	}
	bool operator<(PointID point_id) const { return point < point_id; }
	bool operator>(PointID point_id) const { return point >= point_id; }

	bool valid() const { return point.valid(); }
	void invalidate() { point.invalidate(); }
};

struct FSCoordinate {
	PointID x;
	PointID y;
};

// Orientation and Direction etc.

enum class Direction : uint8_t {
	Up = 0,
	Down = 1,
	Left = 2,
	Right = 3
};
const std::array<Direction, 4> Directions = {{
	Direction::Up,
	Direction::Down,
	Direction::Left,
	Direction::Right
}};

enum class Orientation {
	Horizontal = 0,
	Vertical = 1
};
Orientation operator!(Orientation orientation);

// short for: backward-forward direction
enum class BFDirection {
	Backward = 0,
	Forward = 1
};
BFDirection operator!(BFDirection direction);

// some helper functions for Orientation and Direction
Orientation getOrientation(Direction direction);
Direction getForwardDirection(Orientation orientation);
Direction getBackwardDirection(Orientation orientation);
bool isForward(Direction direction);
bool isBackward(Direction direction);
std::array<Direction, 2> getDirections(Orientation orientation);
BFDirection toBFDirection(Direction direction);

//
// Interval
//

struct Interval
{
	distance_t begin;
	distance_t end;

	Interval()
		: begin(1.),
		  end(0.) {}

	Interval(distance_t begin, distance_t end)
		: begin(begin),
		  end(end) {}

	bool operator<(Interval const& other) const {
		return begin < other.begin || (begin == other.begin && end < other.end);
	}

	bool is_empty() const { return begin > end; }
	bool intersects(Interval const& other) const
	{
		if (is_empty() || other.is_empty()) { return false; }

		return (other.begin >= begin && other.begin <= end) ||
			(other.end >= begin && other.end <= end) ||
			(other.begin <= begin && other.end >= end);
	}
};
using Intervals = std::vector<Interval>;

std::ostream& operator<<(std::ostream& out, const Interval& interval);


// Data Types for FrechetLight:

//NOTES: CPoint will be integral part + rational [0,1] number given as sqrt_extension 

class CPoint {
private:
	PointID point;
	distance_t fraction;

	void normalize() {
		assert(fraction >= 0. && fraction <= 1.);
		if (fraction == 1.) {
			fraction = 0.;
			++point;
		}
	}
public:
	CPoint(PointID point, distance_t fraction)
		: point(point), fraction(fraction)
	{
			normalize();
	}
	CPoint() : point(std::numeric_limits<PointID::IDType>::max()), fraction(0.) {}
	
	PointID getPoint() const { return point; } 
	distance_t getFraction() const { return fraction; } 
	distance_t convert() const { return (distance_t) point + fraction; }
	void setPoint(PointID point) { this->point = point; }
	void setFraction(distance_t frac) { fraction = frac; normalize(); }
		
	bool operator<(CPoint const& other) const {
		return point < other.point || (point == other.point && fraction < other.fraction);
	}
	bool operator<=(CPoint const& other) const {
		return point < other.point || (point == other.point && fraction <= other.fraction);
	}
	bool operator>(CPoint const& other) const {
		return point > other.point || (point == other.point && fraction > other.fraction);
	}
	bool operator>=(CPoint const& other) const {
		return point > other.point || (point == other.point && fraction >= other.fraction);
	}
	bool operator==(CPoint const& other) const {
		return point == other.point && fraction == other.fraction;
	}
	bool operator!=(CPoint const& other) const {
		return point != other.point or fraction != other.fraction;
	}
	bool operator<(PointID other) const {
		return point < other;
	}
	bool operator>(PointID other) const {
		return point > other || (point == other && fraction > 0.);
	}
	bool operator<=(PointID other) const {
		return point < other || (point == other && fraction == 0.);
	}
	bool operator>=(PointID other) const {
		return point >= other;
	}
	bool operator==(PointID other) const {
		return point == other && fraction == 0.;
	}
	bool operator!=(size_t other) const {
		return !(point == other);
	}
	CPoint ceil() const {
		return fraction > 0 ? CPoint(point + 1, 0.) : CPoint(point, 0.);
	}
	CPoint floor() const {
		return CPoint(point, 0.);
	}
	std::string to_string() const { 
	  //return std::to_string( (double) point + fraction); 
	  std::stringstream stream;
	  stream << std::fixed << std::setprecision(10) << (double) point + fraction;
	  return stream.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const CPoint& p);
};

struct CInterval;
using CIntervals = std::vector<CInterval>;
using CIntervalsID = ID<CIntervals>;
using CIntervalID = std::size_t;
using CIntervalIDs = std::vector<CIntervalID>;

using CPoints = std::vector<CPoint>;

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

using CurveID = std::size_t;
using CurveIDs = std::vector<CurveID>;

struct CInterval
{
	CPoint begin;
	CPoint end;

	const CInterval* reach_parent = nullptr; 
	CPoint fixed = CPoint(std::numeric_limits<PointID::IDType>::max(),0.);
	CurveID fixed_curve = -1;

	CPosition getLowerRightPos() const { 
		if (fixed_curve == 0) {
			CPosition ret = {{fixed, begin}}; 
			return ret;
		} else {
			CPosition ret = {{end, fixed}};
			return ret;
		}
	}
	CPosition getUpperLeftPos() const { 
		if (fixed_curve == 0) {
			CPosition ret = {{fixed, end}}; 
			return ret;
		} else {
			CPosition ret = {{begin, fixed}};
			return ret;
		}
	}

	CInterval(CPoint begin, CPoint end, CPoint fixed, CurveID fixed_curve)
		: begin(begin), end(end), fixed(fixed), fixed_curve(fixed_curve) {}

	CInterval()
		: begin(std::numeric_limits<PointID::IDType>::max(), 0.),
		  end(std::numeric_limits<PointID::IDType>::lowest(), 0.) {}

	CInterval(CInterval const& other) = default;

	CInterval(CPoint const& begin, CPoint const& end)
                : begin(begin), end(end) {}

	CInterval(PointID point1, distance_t fraction1, PointID point2, distance_t fraction2)
		: begin(point1, fraction1), end(point2, fraction2) {}

	CInterval(PointID begin, PointID end)
		: begin(begin, 0.), end(end, 0.) {}

	
	bool operator<(CInterval const& other) const {
		return begin < other.begin || (begin == other.begin && end < other.end);
	}

	bool is_empty() const { return end < begin; }
	void make_empty() {
		begin = {std::numeric_limits<PointID::IDType>::max(), 0.};
		end = {std::numeric_limits<PointID::IDType>::lowest(), 0.};
	}
	void clamp(CPoint const& min, CPoint const& max) {
		begin = std::max(min, begin);
		end = std::min(max, end);
	}
};

std::ostream& operator<<(std::ostream& out, const CInterval& interval);

class IntersectionAlgorithm
{
public:
	static constexpr distance_t eps = 1e-8;
	
   /*
    * Returns which section of the line segment from line_start to line_end is inside the circle given by circle_center and radius.
    * If the circle and line segment do not intersect, the result is the empty Interval (and outer is the empty Interval, too).
	* Otherwise the result is an interval [x,y] such that the distance at x and at y is at most the radius, i.e., [x,y] is a subset of the free interval. 
	* The optional output "outer" is an interval strictly containing the free interval.
	* In other words, "outer" is an interval [x',y'] containing [x,y] such that x-x', y'-y <= eps and:
	* If x = 0 then x' = -eps, while if x > 0 then the distance at x' is more than the radius.
	* If y = 1 then y' = 1+eps, while if y < 1 then the distance at y' is more than the radius.
    */
	static Interval intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end, Interval * outer = nullptr);
private:
	IntersectionAlgorithm() {} // Make class static-only
	static inline bool smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr);
	static inline distance_t distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center);

	static constexpr distance_t save_eps = 0.5 * eps;
	static constexpr distance_t save_eps_half = 0.25 * eps;
};

// Ellipse
struct Ellipse
{
	Point center;
	distance_t width;
	distance_t height;
	double alpha;

	void invalidate() { width = -1.; height = -1.; }
	bool is_valid() { return width >= 0 && height >= 0; }
};

Ellipse segmentsToEllipse(Point const& a1, Point const& b1, Point const& a2, Point const& b2, distance_t distance);

// Circle
struct Circle
{
	Point center;
	distance_t radius;

	Circle() = default;
	Circle(Point const& p, distance_t r) : center(p), radius(r) {}
};

namespace
{

template<typename T>
T pow2(T d) { return std::pow(d, 2); }

} // end anonymous namespace

//
// Point
//

Point& Point::operator-=(const Point& point)
{
    x -= point.x;
    y -= point.y;
    return *this;
}

Point Point::operator-(const Point& point) const
{
    auto result = *this;
	result -= point; 

    return result;
}

Point& Point::operator+=(const Point& point)
{
    x += point.x;
    y += point.y;
    return *this;
}

Point Point::operator+(const Point& point) const
{
    auto result = *this;
	result += point; 

    return result;
}

Point Point::operator*(const distance_t mult) const
{
	Point res;
	res.x = mult * this->x;
	res.y = mult * this->y;
    return res;
}

Point& Point::operator/=(distance_t distance)
{
    x /= distance;
    y /= distance;
    return *this;
}

Point Point::operator/(distance_t distance)
{
    return Point(x/distance, y/distance);
}

bool Point::operator==(Point const& other) const
{
	return x == other.x && y == other.y;
}

bool Point::operator!=(Point const& other) const
{
	return !(*this == other);
}

distance_t Point::dist_sqr(const Point& point) const
{
    return pow2(x - point.x) + pow2(y - point.y);
}

distance_t Point::dist(const Point& point) const
{
    return std::sqrt(dist_sqr(point));
}

std::ostream& operator<<(std::ostream& out, const Point& p)
{
    out << std::setprecision (15)
		<< "(" << p.x << ", " << p.y << ")";

    return out;
}

std::ostream& operator<<(std::ostream& out, const CPoint& p)
{
    out << std::setprecision (15)
		<< "(" << (size_t) p.point << " + " << p.fraction << ")";

    return out;
}

//
// Interval
//

std::ostream& operator<<(std::ostream& out, const Interval& interval)
{
    out << std::setprecision (15)
		<< "(" << interval.begin << ", " << interval.end << ")";

    return out;
}

std::ostream& operator<<(std::ostream& out, const CInterval& interval)
{
    out << std::setprecision (15)
		<< "(" << interval.begin << ", " << interval.end << ")";

    return out;
}

//
// intersection_interval
//

inline bool IntersectionAlgorithm::smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr) {
	return circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate) <= radius_sqr;
}

inline distance_t IntersectionAlgorithm::distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center) {
	return circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate);
}

Interval IntersectionAlgorithm::intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end, Interval * outer /* = nullptr*/)
{
    // The line can be represented as line_start + lambda * v
    const Point v = line_end - line_start;
	const distance_t rad_sqr = radius * radius;
	
    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 = radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2) + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) = radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 * line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, 
	// let b := line_start.x * v.x + line_start.y * v.y, 
	// let c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
	
    const distance_t a = pow2(v.x) + pow2(v.y);
    const distance_t b = (line_start.x - circle_center.x) * v.x + (line_start.y - circle_center.y) * v.y;
    const distance_t c = pow2(line_start.x - circle_center.x) + pow2(line_start.y - circle_center.y) - pow2(radius);

	distance_t mid = - b / a;
    distance_t discriminant = pow2(mid) - c / a;

	const bool smallDistAtZero = smallDistanceAt(0., line_start, line_end, circle_center, rad_sqr);
	const bool smallDistAtOne = smallDistanceAt(1., line_start, line_end, circle_center, rad_sqr);
	bool smallDistAtMid = smallDistanceAt(mid, line_start, line_end, circle_center, rad_sqr);
	
	if (smallDistAtZero && smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(-eps, 1. + eps); }
		return Interval(0, 1);
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
		return Interval(); // no intersection;
    }
	
	if (mid <= 0. and !smallDistAtZero) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
	}
	if (mid >= 1. and !smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
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
	
	return Interval{ begin, end };
}

Ellipse segmentsToEllipse(Point const& a1, Point const& b1, Point const& a2, Point const& b2, distance_t distance)
{
	Ellipse e;

	auto const pi = std::atan2(0, -1);

	// Check if segments are parallel
	auto dir1 = b1 - a1;
	auto dir2 = b2 - a2;
	dir1 /= (std::sqrt(pow2(dir1.x)+pow2(dir1.y)));
	dir2 /= (std::sqrt(pow2(dir2.x)+pow2(dir2.y)));
	if (std::abs(dir1.x*dir2.x + dir1.y*dir2.y) >= 0.999) {
		e.invalidate();
		return e;
	}

	// First assign coefficients of x^2, y^2, ...
	// Those are calculated by hand
	auto A = pow2(a1.x) - 2*a1.x*b1.x + pow2(a1.y) - 2*a1.y*b1.y + pow2(b1.x) + pow2(b1.y);
	auto B = -2*a1.x*a2.x + 2*a1.x*b2.x - 2*a1.y*a2.y + 2*a1.y*b2.y + 2*a2.x*b1.x + 2*a2.y*b1.y - 2*b1.x*b2.x - 2*b1.y*b2.y;
	auto C = pow2(a2.x) - 2*a2.x*b2.x + pow2(a2.y) - 2*a2.y*b2.y + pow2(b2.x) + pow2(b2.y);
	auto D = 2*a1.x*b1.x - 2*a1.x*b2.x + 2*a1.y*b1.y - 2*a1.y*b2.y - 2*pow2(b1.x) + 2*b1.x*b2.x - 2*pow2(b1.y) + 2*b1.y*b2.y;
	auto E = -2*a2.x*b1.x + 2*a2.x*b2.x - 2*a2.y*b1.y + 2*a2.y*b2.y + 2*b1.x*b2.x + 2*b1.y*b2.y - 2*pow2(b2.x) - 2*pow2(b2.y);
	auto F = pow2(b1.x) - 2*b1.x*b2.x + pow2(b1.y) - 2*b1.y*b2.y + pow2(b2.x) + pow2(b2.y) - pow2(distance);

	// This should not fail if they are not parallel
	assert(pow2(B) - 4*A*C <= 0.0);

	// Now convert those to the form of the Ellipse type
	// See: https://en.wikipedia.org/wiki/Ellipse#Canonical_form
	e.center.x = (2*C*D - B*E)/(pow2(B) - 4*A*C);
	e.center.y = (2*A*E - B*D)/(pow2(B) - 4*A*C);
	e.width = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C+std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);
	e.height = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C-std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);

	if (B == 0) {
		e.alpha = A < C ? 0 : pi/2;
	}
	else {
		e.alpha = std::atan((C-A-std::sqrt(pow2(A-C) + pow2(B)))/B);
	}
	e.alpha = e.alpha*360/(2*pi);

	return e;
}
