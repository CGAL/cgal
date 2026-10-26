/*
Copyright (C) 2026  TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Point_2.h>

#include <vector>

namespace CGAL {
template <class K> 
struct MultiPoint {
	std::vector<Point_2<K>> points;

	MultiPoint<K> transform(const CGAL::Aff_transformation_2<K>& trans) const {
		MultiPoint<K> transformed;
		for (const auto& p : points) {
			transformed.points.push_back(p.transform(trans));
		}
		return transformed;
	}

	Bbox_2 bbox() const {
		return CGAL::bbox_2(points.begin(), points.end());
	}
};
} // namespace CGAL
