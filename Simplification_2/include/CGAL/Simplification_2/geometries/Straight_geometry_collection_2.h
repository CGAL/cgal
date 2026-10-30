#pragma once

#include "Straight_geometry_2.h"
#include "transform_helpers.h"

#include <CGAL/utility.h>

#include <variant>
#include <vector>

namespace CGAL {
template <class K, class Data, class Container>
class Straight_geometry_collection_2;

namespace internal {
template <class K, class Data, class Container>
class StraightGeometryCollectionHelper;

template <class Geometry> struct Handle {
  public:
	size_t index;
	Handle(size_t index) : index(index) {}

	template <class K, class Data, class Container> friend class CGAL::Straight_geometry_collection_2;
	template <class K, class Data, class Container> friend class StraightGeometryCollectionHelper;
};

template <class K, class Data, class Container> class StraightGeometryCollectionHelper {
  protected:
	std::vector<Straight_geometry_2<K, Container>> m_geometries;
	using Data_container = std::vector<Data>;
	Data_container m_data;

  public:
	Handle<CGAL::Polygon_2<K, Container>> insert(const CGAL::Polygon_2<K, Container>& polygon, Data data) {
		m_geometries.emplace_back(polygon);
		m_data.push_back(std::move(data));
		return Handle<CGAL::Polygon_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<CGAL::Polygon_with_holes_2<K, Container>>
	insert(const CGAL::Polygon_with_holes_2<K, Container>& polygonWithHoles, Data data) {
		m_geometries.emplace_back(polygonWithHoles);
		m_data.push_back(std::move(data));
		return Handle<CGAL::Polygon_with_holes_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<Multipolygon_with_holes_2<K, Container>>
	insert(const Multipolygon_with_holes_2<K, Container>& multipolygon, Data data) {
		m_geometries.emplace_back(multipolygon);
		m_data.push_back(std::move(data));
		return Handle<Multipolygon_with_holes_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<Multipolyline_2<K>> insert(const Multipolyline_2<K>& polyline, Data data) {
		m_geometries.emplace_back(polyline);
		m_data.push_back(std::move(data));
		return Handle<Multipolyline_2<K>>{m_geometries.size() - 1};
	}

	Handle<Polyline_2<K>> insert(const Polyline_2<K>& polyline, Data data) {
		m_geometries.emplace_back(polyline);
		m_data.push_back(std::move(data));
		return Handle<Polyline_2<K>>{m_geometries.size() - 1};
	}

	Handle<Multipoint_2<K>> insert(const Multipoint_2<K>& point, Data data) {
		m_geometries.emplace_back(point);
		m_data.push_back(std::move(data));
		return Handle<Multipoint_2<K>>{m_geometries.size() - 1};
	}

	Handle<CGAL::Point_2<K>> insert(const CGAL::Point_2<K>& point, Data data) {
		m_geometries.emplace_back(point);
		m_data.push_back(std::move(data));
		return Handle<CGAL::Point_2<K>>{m_geometries.size() - 1};
	}

	template <typename T> 
	typename Data_container::reference // for std::vector<bool> this returns a bool by value.
	get_data(Handle<T> handle) {
		return m_data[handle.index];
	}

	template <typename T> 
	typename Data_container::const_reference
	get_data(Handle<T> handle) const {
		return m_data[handle.index];
	}

	Data_container& data() {
		return m_data;
	}

	const Data_container& data() const {
		return m_data;
	}

	Data_container::iterator data_begin() {
		return m_data.begin();
	}

	Data_container::iterator data_end() {
		return m_data.end();
	}

	Data_container::iterator data_begin() const {
		return m_data.begin();
	}

	Data_container::iterator data_end() const {
		return m_data.end();
	}
};

template <class K, class Container>
class StraightGeometryCollectionHelper<K, std::monostate, Container> {
  protected:
	std::vector<Straight_geometry_2<K, Container>> m_geometries;

  public:
	Handle<CGAL::Polygon_2<K, Container>> insert(const CGAL::Polygon_2<K, Container>& polygon) {
		m_geometries.emplace_back(polygon);
		return Handle<CGAL::Polygon_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<CGAL::Polygon_with_holes_2<K, Container>>
	insert(const CGAL::Polygon_with_holes_2<K, Container>& polygonWithHoles) {
		m_geometries.emplace_back(polygonWithHoles);
		return Handle<CGAL::Polygon_with_holes_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<Multipolygon_with_holes_2<K, Container>>
	insert(const Multipolygon_with_holes_2<K, Container>& multipolygon) {
		m_geometries.emplace_back(multipolygon);
		return Handle<Multipolygon_with_holes_2<K, Container>>{m_geometries.size() - 1};
	}

	Handle<Multipolyline_2<K>> insert(const Multipolyline_2<K>& polyline) {
		m_geometries.emplace_back(polyline);
		return Handle<Multipolyline_2<K>>{m_geometries.size() - 1};
	}

	Handle<Polyline_2<K>> insert(const Polyline_2<K>& polyline) {
		m_geometries.emplace_back(polyline);
		return Handle<Polyline_2<K>>{m_geometries.size() - 1};
	}

	Handle<Multipoint_2<K>> insert(const Multipoint_2<K>& point) {
		m_geometries.emplace_back(point);
		return Handle<Multipoint_2<K>>{m_geometries.size() - 1};
	}

	Handle<CGAL::Point_2<K>> insert(const CGAL::Point_2<K>& point) {
		m_geometries.emplace_back(point);
		return Handle<CGAL::Point_2<K>>{m_geometries.size() - 1};
	}
};
} // namespace internal

template <class K, class Data = std::monostate, class Container = std::vector<CGAL::Point_2<K>>>
class Straight_geometry_collection_2
    : public internal::StraightGeometryCollectionHelper<K, Data, Container> {
  public:
	using Geometry_container = std::vector<Straight_geometry_2<K, Container>>;
	template <class Geometry> using Handle = internal::Handle<Geometry>;

	template <typename T> T& get_geometry(Handle<T> handle) {
		return std::get<T>(this->m_geometries.at(handle.index));
	}

	template <typename T> const T& get_geometry(Handle<T> handle) const {
		return std::get<T>(this->m_geometries.at(handle.index));
	}

	Geometry_container& geometries() {
		return this->m_geometries;
	}

	const Geometry_container& geometries() const {
		return this->m_geometries;
	}

	Geometry_container::iterator geometries_begin() {
		return this->m_geometries.begin();
	}

	Geometry_container::iterator geometries_end() {
		return this->m_geometries.end();
	}

	Geometry_container::const_iterator geometries_begin() const {
		return this->m_geometries.begin();
	}

	Geometry_container::const_iterator geometries_end() const {
		return this->m_geometries.end();
	}

	Straight_geometry_collection_2<K, Data, Container> transform(const CGAL::Aff_transformation_2<K>& trans) const {
		Straight_geometry_collection_2<K, Data, Container> transformed(*this);
		using Geometry = Straight_geometry_2<K, Container>;
		for (auto& g : transformed.m_geometries) {
			g = std::visit(
			    overloaded{
			        [&](const CGAL::Polygon_with_holes_2<K, Container>& pwh) -> Geometry { return CGAL::transform(trans, pwh); },
			        [&](const CGAL::Multipolygon_with_holes_2<K, Container>& mp) -> Geometry { return CGAL::transform(trans, mp); },
			        [&](const CGAL::Polygon_2<K, Container>& p) -> Geometry { return CGAL::transform(trans, p); },
			        [&](const auto& someG) -> Geometry { return someG.transform(trans); }},
			    g);
		}
		return transformed;
	}

	CGAL::Bbox_2 bbox() const {
		if (this->m_geometries.empty()) return {};
		CGAL::Bbox_2 box = std::visit([](const auto& geom) { return geom.bbox(); }, *geometries_begin());
		for (auto git = ++geometries_begin(); git != geometries_end(); ++git) {
			box += std::visit([](const auto& geom) { return geom.bbox(); }, *git);
		}
		return box;
	}
};
} // namespace CGAL