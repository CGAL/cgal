// Copyright (c) 2026  Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_AABB_TWO_TREES_TRAVERSAL_TRAITS_H
#define CGAL_AABB_TWO_TREES_TRAVERSAL_TRAITS_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree/internal/Primitive_helper.h>
#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <optional>

namespace CGAL {

namespace internal { namespace AABB_tree {

template<typename AABBTraits1, typename AABBTraits2, typename OutputIterator>
class Two_trees_listing_intersecting_primitives_traits
{
  typedef typename AABBTraits1::Primitive Primitive1;
  typedef typename AABBTraits2::Primitive Primitive2;
  typedef ::CGAL::AABB_node<AABBTraits1> Node1;
  typedef ::CGAL::AABB_node<AABBTraits2> Node2;

  template<bool in_order, typename Value>
  class WrapOutputIterator
  {
    Value first;
    OutputIterator out;
  public:
    WrapOutputIterator(Value first_, OutputIterator out_): first(first_), out(out_){}

    WrapOutputIterator& operator=(Value second){
      if constexpr(in_order)
        out = std::make_pair(first, second);
      else
        out = std::make_pair(second, first);
      return *this;
    }
    WrapOutputIterator& operator*(){ return *this; }
    WrapOutputIterator& operator++(){ ++out; return *this; }
    WrapOutputIterator& operator++(int){ auto tmp = *this; ++out; return tmp;  }
    WrapOutputIterator& operator+(int d){ out += d; return *this; }
  };

public:
  Two_trees_listing_intersecting_primitives_traits(const AABBTraits1& traits1, const AABBTraits2& traits2, OutputIterator out_)
    : m_traits1(traits1), m_traits2(traits2), out(out_)
  {}

  bool go_further() const {
    return true;
  }

  // If true, the next step of traversal continue on A, if false, the next step of traversal is on B
  template<typename Node_A, typename Node_B>
  bool prefer_A_for_next_step(const Node_A& node_a, const Node_B& node_b, const std::size_t&, const std::size_t&) const {
  #if 1
    auto sq_diag = [](const Bbox_3 &b){
      return (b.xmax()-b.xmin())*(b.xmax()-b.xmin()) + (b.ymax()-b.ymin())*(b.ymax()-b.ymin()) + (b.zmax()-b.zmin())*(b.zmax()-b.zmin());
    };
    return sq_diag(node_a.bbox()) > sq_diag(node_b.bbox());
  #else
    return false;
  #endif
  }

  void intersection(const Primitive1& primitive1, const Node2& node2, std::size_t nb_primitives_2)
  {
    using WrapIterator = WrapOutputIterator<true, typename Primitive1::Id>;
    WrapIterator wrap_out(primitive1.id(), out);
    Listing_primitive_traits<AABBTraits2, typename AABBTraits1::Primitive::Datum, WrapIterator> traits(wrap_out, m_traits2);
    node2.traversal( internal::Primitive_helper<AABBTraits1>::get_datum(primitive1, m_traits1), traits, nb_primitives_2);
  }

  void intersection(const Node1& node1, std::size_t nb_primitives_1, const Primitive2& primitive2)
  {
    using WrapIterator= WrapOutputIterator<false, typename Primitive2::Id>;
    WrapIterator wrap_out(primitive2.id(), out);
    Listing_primitive_traits<AABBTraits1, typename AABBTraits2::Primitive::Datum, WrapIterator> traits(wrap_out, m_traits1);
    node1.traversal( internal::Primitive_helper<AABBTraits2>::get_datum(primitive2, m_traits2), traits, nb_primitives_1);
  }

  bool do_intersect(const Node1& node1, const Node2& node2) const
  {
    return do_overlap(node1.bbox(), node2.bbox());
  }

private:
  const AABBTraits1& m_traits1;
  const AABBTraits2& m_traits2;
  OutputIterator out;
};

template<typename AABBTraits, typename OutputIterator>
class Listing_self_intersecting_primitives_traits
{
  typedef typename AABBTraits::Primitive Primitive;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

  template<bool in_order, typename Value>
  class WrapOutputIterator
  {
    Value first;
    OutputIterator out;
  public:
    WrapOutputIterator(Value first_, OutputIterator out_): first(first_), out(out_){}

    WrapOutputIterator& operator=(Value second){
      if constexpr(in_order)
        out = std::make_pair(first, second);
      else
        out = std::make_pair(second, first);
      return *this;
    }
    WrapOutputIterator& operator*(){ return *this; }
    WrapOutputIterator& operator++(){ ++out; return *this; }
    WrapOutputIterator& operator++(int){ auto tmp = *this; ++out; return *this;  }
    WrapOutputIterator& operator+(int d){ out += d; return *this; }
  };

public:
  Listing_self_intersecting_primitives_traits(const AABBTraits& traits, OutputIterator out_)
    : m_traits(traits), out(out_)
  {}

  bool go_further() const {
    return true;
  }

  // If true, the next step of traversal continue on A, if false, the next step of traversal is on B
  template<typename Node_A, typename Node_B>
  bool prefer_A_for_next_step(const Node_A& node_a, const Node_B& node_b, const std::size_t&, const std::size_t&) const {
    auto sq_diag = [](const Bbox_3 &b){
      return (b.xmax()-b.xmin())*(b.xmax()-b.xmin()) + (b.ymax()-b.ymin())*(b.ymax()-b.ymin()) + (b.zmax()-b.zmin())*(b.zmax()-b.zmin());
    };
    return sq_diag(node_a.bbox()) > sq_diag(node_b.bbox());
  }

  void intersection(const Primitive& primitive1, const Node& node2, std::size_t nb_primitives_2)
  {
    // TODO Since we are in a symetrical case, maybe we can ignore this call
    using WrapIterator = WrapOutputIterator<true, typename Primitive::Id>;
    WrapIterator wrap_out(primitive1.id(), out);
    Listing_distinct_primitive_traits<AABBTraits, WrapIterator> traits(wrap_out, m_traits);
    node2.traversal( primitive1, traits, nb_primitives_2);
  }

  void intersection(const Node& node1, std::size_t nb_primitives_1, const Primitive& primitive2)
  {
    using WrapIterator= WrapOutputIterator<false, typename Primitive::Id>;
    WrapIterator wrap_out(primitive2.id(), out);
    Listing_distinct_primitive_traits<AABBTraits, WrapIterator> traits(wrap_out, m_traits);
    node1.traversal( primitive2, traits, nb_primitives_1);
  }

  bool do_intersect(const Node& node1, const Node& node2) const
  {
    return do_overlap(node1.bbox(), node2.bbox());
  }

private:
  const AABBTraits& m_traits;
  OutputIterator out;
};

template<typename AABBTraits1, typename AABBTraits2>
class Two_trees_do_intersect_traits
{
  typedef typename AABBTraits1::Primitive Primitive1;
  typedef typename AABBTraits2::Primitive Primitive2;
  typedef ::CGAL::AABB_node<AABBTraits1> Node1;
  typedef ::CGAL::AABB_node<AABBTraits2> Node2;

public:
  Two_trees_do_intersect_traits(const AABBTraits1& traits1, const AABBTraits2& traits2)
    : m_traits1(traits1), m_traits2(traits2), m_is_found(false)
  {}

  bool go_further() const {
    return !m_is_found;
  }

  // If true, the next step of traversal continue on A, if false, the next step of traversal is on B
  template<typename Node_A, typename Node_B>
  bool prefer_A_for_next_step(const Node_A& node_a, const Node_B& node_b, const std::size_t&, const std::size_t&) const {
  #if 1
    auto sq_diag = [](const Bbox_3 &b){
      return (b.xmax()-b.xmin())*(b.xmax()-b.xmin()) + (b.ymax()-b.ymin())*(b.ymax()-b.ymin()) + (b.zmax()-b.zmin())*(b.zmax()-b.zmin());
    };
    return sq_diag(node_a.bbox()) > sq_diag(node_b.bbox());
  #else
    return false;
  #endif
  }

  void intersection(const Primitive1& primitive1, const Node2& node2, std::size_t nb_primitives_2)
  {
    Do_intersect_traits<AABBTraits2, typename AABBTraits1::Primitive::Datum> traits(m_traits2);
    node2.traversal( internal::Primitive_helper<AABBTraits1>::get_datum(primitive1, m_traits1), traits, nb_primitives_2);
    if(traits.is_intersection_found())
      m_is_found = true;
  }

  void intersection(const Node1& node1, std::size_t nb_primitives_1, const Primitive2& primitive2)
  {
    Do_intersect_traits<AABBTraits1, typename AABBTraits2::Primitive::Datum> traits(m_traits1);
    node1.traversal( internal::Primitive_helper<AABBTraits2>::get_datum(primitive2, m_traits2), traits, nb_primitives_1);
    if(traits.is_intersection_found())
      m_is_found = true;
  }

  bool do_intersect(const Node1& node1, const Node2& node2) const
  {
    return do_overlap(node1.bbox(), node2.bbox());
  }

  bool is_intersection_found() const { return m_is_found; }

private:
  const AABBTraits1& m_traits1;
  const AABBTraits2& m_traits2;
  bool m_is_found;
};

namespace experimental{

template<typename AABBTraits1, typename AABBTraits2, typename OutputIterator>
class Two_trees_intersecting_nodes_traits
{
  typedef typename AABBTraits1::Primitive Primitive1;
  typedef typename AABBTraits2::Primitive Primitive2;
  typedef ::CGAL::AABB_node<AABBTraits1> Node1;
  typedef ::CGAL::AABB_node<AABBTraits2> Node2;

public:
  Two_trees_intersecting_nodes_traits(const AABBTraits1& traits1, const AABBTraits2& traits2, OutputIterator out_)
    : m_traits1(traits1), m_traits2(traits2), out(out_)
  {}

  bool go_further() const {
    return true;
  }

  void intersection(const Node1& node1, const Node2& node2)
  {
    *out++ = std::make_pair(&node1, &node2);
  }

  bool do_intersect(const Node1& node1, const Node2& node2) const
  {
    return do_overlap(node1.bbox(), node2.bbox());
  }

private:
  const AABBTraits1& m_traits1;
  const AABBTraits2& m_traits2;
  OutputIterator out;
};

}


}}} // end namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
