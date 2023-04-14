#pragma once

#include "defs.h"
#include "geometry_basics.h"
#include "id.h"

// TODO: should probably make this a more general type in CGAL
// Represents a trajectory. Additionally to the points given in the input file,
// we also store the length of any prefix of the trajectory.
class Curve
{
public:
    Curve() = default;
    Curve(const Points& points);

    std::size_t size() const { return points.size(); }
	bool empty() const { return points.empty(); }
    Point const& operator[](PointID i) const { return points[i]; }
    bool operator==(Curve const& other) const {
		return std::equal(points.cbegin(), points.cend(), other.points.cbegin(), other.points.cend());
	}
    bool operator!=(Curve const& other) const {
		return !(*this == other);
	}
	Point interpolate_at(CPoint const& pt) const  {
		assert(pt.getFraction() >= 0. && pt.getFraction() <= 1.);
		assert((pt.getPoint() < points.size()-1 || (pt.getPoint() == points.size()-1 && pt.getFraction() == 0.)));
		return pt.getFraction() == 0. ? points[pt.getPoint()] : points[pt.getPoint()]*(1.-pt.getFraction()) + points[pt.getPoint()+1]*pt.getFraction();
	}
	Point interpolate_at(PointID const& pt) const  { return points[pt]; }
    distance_t curve_length(PointID i, PointID j) const
		{ return prefix_length[j] - prefix_length[i]; }

    Point front() const { return points.front(); }
    Point back() const { return points.back(); }

    void push_back(Point const& point);

	Points::const_iterator begin() { return points.begin(); }
	Points::const_iterator end() { return points.end(); }
	Points::const_iterator begin() const { return points.cbegin(); }
	Points::const_iterator end() const { return points.cend(); }
	
	std::string filename;

	struct ExtremePoints { distance_t min_x, min_y, max_x, max_y; };
	ExtremePoints const& getExtremePoints() const;
	distance_t getUpperBoundDistance(Curve const& other) const;

private:
    Points points;
    std::vector<distance_t> prefix_length;
	ExtremePoints extreme_points = {
		std::numeric_limits<distance_t>::max(), std::numeric_limits<distance_t>::max(),
		std::numeric_limits<distance_t>::lowest(), std::numeric_limits<distance_t>::lowest()
	};
};
using Curves = std::vector<Curve>;

std::ostream& operator<<(std::ostream& out, const Curve& curve);

Curve::Curve(const Points& points)
	: points(points), prefix_length(points.size())
{
	if (points.empty()) { return; }

	auto const& front = points.front();
	extreme_points = { front.x, front.y, front.x, front.y };
	prefix_length[0] = 0;

	for (PointID i = 1; i < points.size(); ++i)
	{
		auto segment_distance = points[i - 1].dist(points[i]);
		prefix_length[i] = prefix_length[i - 1] + segment_distance;

		extreme_points.min_x = std::min(extreme_points.min_x, points[i].x);
		extreme_points.min_y = std::min(extreme_points.min_y, points[i].y);
		extreme_points.max_x = std::max(extreme_points.max_x, points[i].x);
		extreme_points.max_y = std::max(extreme_points.max_y, points[i].y);
	}
}

void Curve::push_back(Point const& point)
{
	if (prefix_length.size()) {
		auto segment_distance = points.back().dist(point);
		prefix_length.push_back(prefix_length.back() + segment_distance);
	}
	else {
		prefix_length.push_back(0);
	}

	extreme_points.min_x = std::min(extreme_points.min_x, point.x);
	extreme_points.min_y = std::min(extreme_points.min_y, point.y);
	extreme_points.max_x = std::max(extreme_points.max_x, point.x);
	extreme_points.max_y = std::max(extreme_points.max_y, point.y);

	points.push_back(point);
}

auto Curve::getExtremePoints() const -> ExtremePoints const&
{
	return extreme_points;
}

distance_t Curve::getUpperBoundDistance(Curve const& other) const
{
	auto const& extreme1 = this->getExtremePoints();
	auto const& extreme2 = other.getExtremePoints();

	Point min_point{ std::min(extreme1.min_x, extreme2.min_x),
		std::min(extreme1.min_y, extreme2.min_y) };
	Point max_point = { std::max(extreme1.max_x, extreme2.max_x),
		std::max(extreme1.max_y, extreme2.max_y) };

	return min_point.dist(max_point);
}

std::ostream& operator<<(std::ostream& out, const Curve& curve)
{
    out << "[";
	for (auto const& point: curve) {
		out << point << ", ";
	}
    out << "]";

    return out;
}
