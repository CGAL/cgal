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

#include <concepts>

namespace CGAL::internal {

template <class T> concept GraphTraits_2 = requires {
	T::historic;
	requires std::is_same_v<bool const, decltype(T::historic)>;

	T::oriented;
	requires std::is_same_v<bool const, decltype(T::oriented)>;

	T::sorted;
	requires std::is_same_v<bool const, decltype(T::sorted)>;

	T::decomposed;
	requires std::is_same_v<bool const, decltype(T::decomposed)>;
};

template <class T> concept HasPathData = requires {
	typename T::path_data;
};

template<bool H, bool O, bool S, bool D, typename PD> struct CustomGraphTraits {
	static constexpr bool historic = H;
	static constexpr bool oriented = O;
	static constexpr bool sorted = S;
	static constexpr bool decomposed = D;
	using PathData = PD;
};
static_assert(GraphTraits_2<CustomGraphTraits<true, true, true, true, std::monostate>>);

template <bool H> struct PlainGraph {
	static constexpr bool historic = H;
	static constexpr bool oriented = false;
	static constexpr bool sorted = false;
	static constexpr bool decomposed = false;
	using PathData = std::monostate;
};
static_assert(GraphTraits_2<PlainGraph<true>>);

template<bool H> 
struct OrientedGraph {
	static constexpr bool historic = H;
	static constexpr bool oriented = true;
	static constexpr bool sorted = false;
	static constexpr bool decomposed = false;
	using PathData = std::monostate;
};
static_assert(GraphTraits_2<OrientedGraph<true>>);

template <bool H, typename PD> struct DecomposedGraph {
	static constexpr bool historic = H;
	static constexpr bool oriented = true;
	static constexpr bool sorted = false;
	static constexpr bool decomposed = true;
	using PathData= PD;
};
static_assert(GraphTraits_2<DecomposedGraph<true,std::monostate>>);

} // namespace CGAL::internal
