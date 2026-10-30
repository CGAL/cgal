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
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <iterator>
#include <vector>

namespace CGAL {
template <class K>
Polygon_with_holes_2<K> transform(const CGAL::Aff_transformation_2<K> &t, const Polygon_with_holes_2<K> &pwh) {
	Polygon_2<K> outerT;
	if (!pwh.is_unbounded()) {
		outerT = transform(t, pwh.outer_boundary());
	}
	std::vector <Polygon_2<K>> holesT;
	for (const auto &h: pwh.holes()) {
		holesT.push_back(transform(t, h));
	}
	return {outerT, holesT.begin(), holesT.end()};
}

template <class K>
Polygon_set_2<K> transform(const CGAL::Aff_transformation_2<K> &t, const Polygon_set_2<K> &ps) {
	Polygon_set_2<K> transformed;
	std::vector<Polygon_with_holes_2<K>> pwhs;
	ps.polygons_with_holes(std::back_inserter(pwhs));
	for (const auto& pwh : pwhs) {
		transformed.insert(transform(t, pwh));
	}
	return transformed;
}
}
