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
