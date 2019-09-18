#pragma once

#include "geometry_basics.h"
#include "curve.h"

class Certificate
{
public:
	Certificate() = default;

	bool isYes() const { 
		assert(isValid());
		return lessThan;
	}

	bool isValid() const { return valid; }

	bool check() const;
	const CPositions& getTraversal() const { return traversal; }

	void addPoint(const CPosition& pos) { traversal.push_back(pos); }
	void setAnswer(bool answer) { lessThan = answer; }
	void setCurves(const Curve * curve1, const Curve * curve2) { curve_pair[0] = curve1; curve_pair[1] = curve2; }
	void setDistance(distance_t distance) { dist = distance; dist_sqr = std::pow(distance, 2); }
	void validate() { valid = true; };
	void reset() { valid = false; traversal.clear(); }
	void clear() { traversal.clear(); }

	void dump_certificate() const;

private:
	//if YES certificate: (traversal1[0], traversal2[0]), ..., (traversal1[T], traversal2[T]) should be feasible traversal of curve1 and curve2 (interpolate between discrete points) 
	//if NO certificate: certificate should be as described in certificate-outline.pdf

	CPositions traversal;
	std::array<const Curve *, 2> curve_pair;
	distance_t dist, dist_sqr;

	bool lessThan;
	bool valid = false;

	bool feasible(const CPosition & pt) const;
	bool feasible(const CPoint &pt1, const CPoint &pt2) const;
	bool nonEmpty(CurveID fixed_curve, const CPoint& fixed_point, const CPoint& start_pt, const CPoint& end_point) const;
};
