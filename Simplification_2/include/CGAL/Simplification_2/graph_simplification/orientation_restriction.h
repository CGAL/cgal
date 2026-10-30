#pragma once

#include <CGAL/Simplification_2/core.h>

#include <functional>
#include <initializer_list>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace CGAL::internal {

	// Restricts the orientations of all edges of the graph to the given orientations, by replacing edges with other orientations by staircases.
	// An orientation includes its opposite direction. Edges longer than lambda times the longest edge are first subdivided;
	// eps is the tolerance used for the staircases. The progress callback receives a step label, a counter and a maximum.

	// The orientations are unit vectors with angles in [0, pi), sorted by angle.
	template<class Graph> requires Graph::Graph_traits::sorted && Graph::Graph_traits::oriented
	void restrict_orientations(Graph& graph, std::vector<Vector<Inexact>> orientations, Number<Inexact> lambda = 1, Number<Inexact> eps = 0.1, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt);

	// The orientations are count evenly spaced angles, starting at initial_angle (in radians).
	template<class Graph> requires Graph::Graph_traits::sorted && Graph::Graph_traits::oriented
	void restrict_orientations(Graph& graph, int count, Number<Inexact> initial_angle = 0, Number<Inexact> lambda = 1, Number<Inexact> eps = 0.1, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt);

	// The orientations are angles in radians, in [0, pi) and sorted.
	template<class Graph> requires Graph::Graph_traits::sorted && Graph::Graph_traits::oriented
	void restrict_orientations(Graph& graph, std::vector<Number<Inexact>> angles, Number<Inexact> lambda = 1, Number<Inexact> eps = 0.1, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt);

	// The orientations are angles in radians, in [0, pi) and sorted.
	template<class Graph> requires Graph::Graph_traits::sorted && Graph::Graph_traits::oriented
	void restrict_orientations(Graph& graph, std::initializer_list<Number<Inexact>> angles, Number<Inexact> lambda = 1, Number<Inexact> eps = 0.1, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt);

}

#include "orientation_restriction.hpp"
