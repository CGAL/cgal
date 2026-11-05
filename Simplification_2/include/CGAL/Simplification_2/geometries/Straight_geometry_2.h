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

#include "Polyline_2.h"
#include "Multipolyline_2.h"
#include "Multipoint_2.h"
#include <CGAL/Multipolygon_with_holes_2.h>

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <variant>
#include <vector>

namespace CGAL {
/// Type that serves as an internal representation of the straight (i.e. linear) features of the OGC Simple Feature Access.
/// Its alternatives correspond to the OGC types MultiPolygon, Polygon, LinearRing, MultiLineString, LineString, Point and MultiPoint.
template <class K, class Container = std::vector<Point_2<K>>>
using Straight_geometry_2 = std::variant<Multipolygon_with_holes_2<K, Container>, Polygon_with_holes_2<K, Container>, Polygon_2<K, Container>,
                                         Multipolyline_2<K>, Polyline_2<K>, Point_2<K>, Multipoint_2<K>>;
}