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

#ifndef CGAL_POLYLINE_D_H
#define CGAL_POLYLINE_D_H

#include <CGAL/basic.h>

// TODO: what is necessary here?
#include <iterator>

namespace CGAL {

template <class Traits,
          class Container = std::vector<typename Traits::Point_d> >
class Polyline_d {
public:
	//
	// typedefs
	//

    using Point_d = typename Traits::Point_d;
    using Segment_d = typename Traits::Segment_d;
    using Bbox_d = typename Traits::Bbox_d;

	using Vertex_const_iterator = typename Container::const_iterator;

	//
	// constructors
	//

	// TODO: make sure to use right constructors in Fréchet code
	Polyline_d() : traits() {}

	Polyline_d(const Traits& traits) : traits(traits) {}

	Polyline_d(const Polyline_d<Traits, Container>& polyline)
	: points(polyline.points), traits(polyline.traits) {
		// TODO: update bbox
	}

	template <class InputIterator>
	Polyline_d(InputIterator first, InputIterator last,
	           Traits p_traits = Traits())
	: d_container(), traits(p_traits) {
		std::copy(first, last, std::back_inserter(d_container));
		// TODO: update bbox
	}

	//
	// member functions
	//

    std::size_t size() const {
		return points.size();
	}

	// TODO: replace renaming of this function
	bool is_empty() const {
		return points.empty();
	}

	const Point_d& vertex(std::size_t i) const {
		CGAL_precondition(i < points.size());
		return *(std::next(points.begin(), i));
	}

	const Point_d& operator[](std::size_t i) const {
		return vertex(i);
	}

	// TODO: need this?
    bool operator==(Polyline_d<Traits, Container> const& other) const {
		return std::equal(points.cbegin(), points.cend(), other.points.cbegin(), other.points.cend());
	}

	// TODO: need this?
    bool operator!=(Polyline_d<Traits, Container> const& other) const {
		return !(*this == other);
	}

	////////////////////
	// // TODO: should probably not be contained in this class and instead moved to Fréchet as these are intersections of circle and segment and therefore squareroot extensions. Seems too specialized for Polyline_d.
    // 
	// Point_d interpolate_at(CPoint const& pt) const {
	//     CGAL_precondition(pt.getFraction() >= 0. && pt.getFraction() <= 1.);
	//     CGAL_precondition((pt.getPoint() < points.size()-1 || (pt.getPoint() == points.size()-1 && pt.getFraction() == 0.)));
    // 
	//     return pt.getFraction() == 0. ? points[pt.getPoint()] : points[pt.getPoint()]*(1.-pt.getFraction()) + points[pt.getPoint()+1]*pt.getFraction();
	// }
    // 
	// // TODO: fix using CGAL stuff
	// Point_d interpolate_at(PointID const& pt) const  { return points[pt]; }
    // 
	// // TODO: replace distance_t
    // distance_t curve_length(PointID i, PointID j) const {
	//     return prefix_length[j] - prefix_length[i];
	// }
    // 
	// // TODO: move to Fréchet code
	// distance_t getUpperBoundDistance(Curve const& other) const;
    // 
	// // TODO: have to handle "filename" member differently in the Fréchet code
	// // -> or don't maintain at all as there should be another identifier
	////////////////////

    Point_d front() const {
		return points.front();
	}
    Point_d back() const {
		return points.back();
	}

	// TODO: make sure to update bbox here
	// TODO: add other functions here!
    void push_back(Point_d const& point);

	Vertex_const_iterator begin() { return points.begin(); }
	Vertex_const_iterator end() { return points.end(); }
	Vertex_const_iterator begin() const { return points.cbegin(); }
	Vertex_const_iterator end() const { return points.cend(); }

	// TODO: for bounding box, either
	// * completely recompute from scratch
	// * maintain bounding box with "valid" flag and update directly on insert, on remove only flag not "valid"
	// * do not include in this class and instead maintain in utility class

	// TODO: use CGAL bounding box here
	// struct ExtremePoints { distance_t min_x, min_y, max_x, max_y; };
	// ExtremePoints const& getExtremePoints() const;
	// TODO: replace renaming of this function
	BBox_d const& bbox() const {
		return bbox;
	}

private:
	Container points;
	Traits traits;
	// TODO: remember to update the bounding box in all member functions.
	BBox_d bbox;
};

} // end of namespace CGAL

#endif // CGAL_POLYLINE_D_H
