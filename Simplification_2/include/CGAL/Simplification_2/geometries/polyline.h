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
#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/squared_distance_2.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

namespace CGAL {
template <class K> 
CGAL::Point_2<K> projection(const CGAL::Segment_2<K>& seg, const CGAL::Point_2<K>& p) {
	Line_2<K> l = seg.supporting_line();
	auto q = l.projection(p);
	auto s = (q - seg.source()) * (seg.target() - seg.source()) / seg.squared_length();
	if (s < 0)
		return seg.source();
	if (s > 1)
		return seg.target();
	return q;
}

template <class Segment_>
class Polyline_Segment_ptr
{
  public:
	typedef Segment_ Segment;
	Polyline_Segment_ptr(Segment const &seg) :m_seg(seg){}
	Segment* operator->() {return &m_seg;}
  private:
	Segment m_seg;
};

template <class K, class InputIterator> class SegmentIterator {
  public:
	using iterator_category = std::input_iterator_tag;
	using value_type = CGAL::Segment_2<K>;
	using difference_type = typename InputIterator::difference_type;
	using pointer = CGAL::Segment_2<K>*;
	using reference = CGAL::Segment_2<K>&;
  private:
	InputIterator m_first_vertex;

  public:
	typedef SegmentIterator<K, InputIterator> Self;

	SegmentIterator(InputIterator first_vertex): m_first_vertex(first_vertex) {};

	CGAL::Segment_2<K> operator*() const {
		auto second_vertex = m_first_vertex;
		++second_vertex;
		return { *m_first_vertex, *second_vertex };
	}

	Polyline_Segment_ptr<CGAL::Segment_2<K>> operator->() const
	{return Polyline_Segment_ptr<CGAL::Segment_2<K>>(operator*());}

	Self& operator++() {
		++m_first_vertex;
		return *this;
	}

	Self operator++(int) {
		Self tmp = *this;
		++*this;
		return tmp;
	}

	bool operator==(const Self& other) const {
		return m_first_vertex == other.m_first_vertex;
	}

	bool operator!=(const Self& other) const {
		return m_first_vertex != other.m_first_vertex;
	}
};

template <class K> class Polyline_2 {
  public:
	typedef typename std::vector<CGAL::Point_2<K>>::const_iterator Vertex_iterator;
	typedef typename std::vector<CGAL::Point_2<K>>::const_reverse_iterator Vertex_reverse_iterator;
	typedef SegmentIterator<K, typename std::vector<typename CGAL::Point_2<K>>::const_iterator> Edge_iterator;
	Polyline_2() = default;

	// Container typedefs and functions so that a polyline is a range of points: e.g. std::back_inserter can be used
	// on a polyline, and it can be read and written as a WKT linestring.
	using value_type = Point_2<K>;
	using reference = Point_2<K>&;
	using const_reference = const Point_2<K>&;
	using iterator = Vertex_iterator;
	using const_iterator = Vertex_iterator;
	using reverse_iterator = Vertex_reverse_iterator;
	using const_reverse_iterator = Vertex_reverse_iterator;
	using difference_type = int;
	using size_type = size_t;

	[[nodiscard]] const_iterator begin() const { return m_points.begin(); }
	[[nodiscard]] const_iterator end() const { return m_points.end(); }
	[[nodiscard]] const_reverse_iterator rbegin() const { return m_points.rbegin(); }
	[[nodiscard]] const_reverse_iterator rend() const { return m_points.rend(); }

	void clear() {
		m_points.clear();
	}

	void resize(size_type n) {
		m_points.resize(n);
	}

	template <class InputIterator> Polyline_2(InputIterator begin, InputIterator end) {
		if (begin == end) {
			throw std::runtime_error("Polyline cannot be empty.");
		}
		std::copy(begin, end, std::back_inserter(m_points));
	}

	Polyline_2(Segment_2<K> segment) : m_points({segment.source(), segment.target()}) {};

	explicit Polyline_2(std::vector<CGAL::Point_2<K>> points): m_points(points) {};

	void push_back(const CGAL::Point_2<K>& p) {
		m_points.push_back(p);
	}

	void insert(Vertex_iterator i, const CGAL::Point_2<K>& p) {
		m_points.insert(i, p);
	}

	void pop_back() {
		m_points.pop_back();
	}

	[[nodiscard]] Vertex_iterator vertices_begin() const { return m_points.begin(); }
	[[nodiscard]] Vertex_iterator vertices_end() const { return m_points.end(); }
	[[nodiscard]] Vertex_reverse_iterator vertices_rbegin() const { return m_points.rbegin(); }
	[[nodiscard]] Vertex_reverse_iterator vertices_rend() const { return m_points.rend(); }
	[[nodiscard]] Edge_iterator edges_begin() const { return { m_points.begin() }; }
	[[nodiscard]] Edge_iterator edges_end() const { return { --m_points.end() }; }
	[[nodiscard]] int num_vertices() const { return m_points.size(); }
	[[nodiscard]] int num_edges() const { return m_points.size() - 1; }
	[[nodiscard]] int size() const { return num_vertices(); }
	[[nodiscard]] CGAL::Point_2<K> vertex(int i) const { return m_points[i]; }
	[[nodiscard]] CGAL::Segment_2<K> edge(int i) const { return *(std::next(edges_begin(), i)); }
	[[nodiscard]] CGAL::Bbox_2 bbox() const { return CGAL::bbox_2(m_points.begin(), m_points.end()); }
	[[nodiscard]] CGAL::Point_2<K> source() const { return m_points.front(); }
	[[nodiscard]] CGAL::Point_2<K> target() const { return m_points.back(); }
	[[nodiscard]] Polyline_2<K> transform(const CGAL::Aff_transformation_2<K>& trans) const {
		Polyline_2<K> transformed;
		for (const auto& pt : m_points) {
			transformed.push_back(pt.transform(trans));
		}
		return transformed;
	}
	CGAL::Point_2<K> nearest(const CGAL::Point_2<K>& p) const {
		CGAL::Point_2<K> closest;
		typename K::FT minDist2 = std::numeric_limits<typename K::FT>::infinity();

		for (auto eit = edges_begin(); eit != edges_end(); ++eit) {
			auto proj = projection(*eit, p);
			auto dist2 = CGAL::squared_distance(p, proj);
			if (dist2 < minDist2) {
				minDist2 = dist2;
				closest = proj;
			}
		}

		return closest;
	}
  private:
	std::vector<CGAL::Point_2<K>> m_points;
};
}
