// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_TYPE_DEFINITIONS_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_TYPE_DEFINITIONS_H

#include <CGAL/license/Mesh_smoothing_3.h>


#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_primitive.h>

namespace CGAL {

namespace Mesh_smoothing_3_internal {

template <class C3t3>
typename C3t3::Triangulation::Geom_traits::Point_3
get_point(typename C3t3::Vertex_handle vh)
{
  using T3=typename C3t3::Triangulation;
  if constexpr (std::is_same_v<typename T3::Geom_traits::Weighted_point_3,
                               typename T3::Vertex::Point>)
    return vh->point().point();
  else
    return vh->point();
}

template<typename C3t3>
struct Facet_to_triangle_property_map
{
    using key_type = typename C3t3::Facet;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Triangle_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Facet_to_triangle_property_map<C3t3>::reference
get(const Facet_to_triangle_property_map<C3t3>&, const typename C3t3::Facet& f)
{
    using Tr = typename C3t3::Triangulation;
    using K = typename Tr::Geom_traits;
    const typename Tr::Cell_handle cell = f.first;
    const int i = f.second;
    const typename K::Point_3& p0 = get_point<C3t3>(cell->vertex((i + 1) % 4));
    const typename K::Point_3& p1 = get_point<C3t3>(cell->vertex((i + 2) % 4));
    const typename K::Point_3& p2 = get_point<C3t3>(cell->vertex((i + 3) % 4));
    return typename K::Triangle_3(p0, p1, p2);
}

template<typename C3t3>
struct Facet_to_point_property_map
{
    using key_type = typename C3t3::Facet;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Point_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Facet_to_point_property_map<C3t3>::reference
get(const Facet_to_point_property_map<C3t3>&, const typename C3t3::Facet& f)
{
    using Tr = typename C3t3::Triangulation;
    const typename Tr::Cell_handle cell = f.first;
    const int i = f.second;
    return get_point<C3t3>(cell->vertex((i + 1) % 4));
}

template<typename C3t3>
struct Edge_to_segment_property_map
{
    using key_type = typename C3t3::Edge;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Segment_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Edge_to_segment_property_map<C3t3>::reference
get(const Edge_to_segment_property_map<C3t3>&, const typename C3t3::Edge& e)
{
    using Tr = typename C3t3::Triangulation;
    using K = typename Tr::Geom_traits;
    const typename Tr::Cell_handle cell = e.first;
    return typename K::Segment_3(get_point<C3t3>(cell->vertex(e.second)),
                                 get_point<C3t3>(cell->vertex(e.third)));
}

template<typename C3t3>
struct Edge_to_point_property_map
{
    using key_type = typename C3t3::Edge;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Point_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Edge_to_point_property_map<C3t3>::reference
get(const Edge_to_point_property_map<C3t3>&, const typename C3t3::Edge& e)
{
    using Tr = typename C3t3::Triangulation;
    const typename Tr::Cell_handle cell = e.first;
    return get_point<C3t3>(cell->vertex(e.second));
}

template <typename C3t3>
using Facet_primitive = CGAL::AABB_primitive<typename C3t3::Facet,
                                                Facet_to_triangle_property_map<C3t3>,
                                                Facet_to_point_property_map<C3t3>,
                                                CGAL::Tag_false,
                                                CGAL::Tag_false>;

template <typename C3t3>
using Facet_traits = CGAL::AABB_traits_3<typename C3t3::Triangulation::Geom_traits, Facet_primitive<C3t3>>;

template <typename C3t3>
using Facet_tree = CGAL::AABB_tree<Facet_traits<C3t3>>;


template <typename C3t3>
using Edge_primitive = CGAL::AABB_primitive<typename C3t3::Edge,
                                                Edge_to_segment_property_map<C3t3>,
                                                Edge_to_point_property_map<C3t3>,
                                                CGAL::Tag_false,
                                                CGAL::Tag_false>;
template <typename C3t3>
using Edge_traits = CGAL::AABB_traits_3<typename C3t3::Triangulation::Geom_traits, Edge_primitive<C3t3>>;

template <typename C3t3>
using Edge_tree = CGAL::AABB_tree<Edge_traits<C3t3>>;



} } // end of CGAL::Mesh_smoothing_3_internal namespace

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_TYPE_DEFINITIONS_H
