#pragma once

#include "geometry_basics.h"

struct Circle
{
	Point center;
	distance_t radius;

	Circle() = default;
	Circle(Point const& point, distance_t radius) : center(point), radius(radius) {}
};

struct LineArc
{
	Point start;
	Point end;

	LineArc() = default;
	LineArc(Point const& start, Point const& end) : start(start), end(end) {}
};

using Circular_arc_point_2 = Point;
