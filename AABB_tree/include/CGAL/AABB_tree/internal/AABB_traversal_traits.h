// Copyright (c) 2008-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_TRAVERSAL_TRAITS_H
#define CGAL_AABB_TRAVERSAL_TRAITS_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <boost/optional.hpp>

namespace CGAL {

namespace internal { namespace AABB_tree {

template <class Value_type, typename Integral_type>
class Counting_output_iterator {
  typedef Counting_output_iterator<Value_type,Integral_type> Self;
  Integral_type* i;
public:
  Counting_output_iterator(Integral_type* i_) : i(i_) {};

  struct Proxy {
    Proxy& operator=(const Value_type&) { return *this; };
  };

  Proxy operator*() {
    return Proxy();
  }

  Self& operator++() {
    ++*i;
    return *this;
  }

  Self& operator++(int) {
    ++*i;
    return *this;
  }
};

//-------------------------------------------------------
// Traits classes for traversal computation
//-------------------------------------------------------
/**
 * @class First_intersection_traits
 */
template<typename AABBTraits, typename Query>
class First_intersection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  typedef
  boost::optional< typename AABBTraits::template Intersection_and_primitive_id<Query>::Type >
  Result;
public:
  First_intersection_traits(const AABBTraits& traits)
    : m_result(), m_traits(traits)
  {}

  bool go_further() const {
    return !m_result;
  }

  void intersection(const Query& query, const Primitive& primitive)
  {
    m_result = m_traits.intersection_object()(query, primitive);
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

  Result result() const { return m_result; }
  bool is_intersection_found() const {
    return m_result;
  }

private:
  Result m_result;
  const AABBTraits& m_traits;
};


/**
 * @class Listing_intersection_traits
 */
template<typename AABBTraits, typename Query, typename Output_iterator>
class Listing_intersection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  Listing_intersection_traits(Output_iterator out_it, const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {}

  constexpr bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    boost::optional< typename AABBTraits::template Intersection_and_primitive_id<Query>::Type >
    intersection = m_traits.intersection_object()(query, primitive);

    if(intersection)
    {
      *m_out_it++ = *intersection;
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;
  const AABBTraits& m_traits;
};


/**
 * @class Listing_primitive_traits
 */
template<typename AABBTraits, typename Query, typename Output_iterator>
class Listing_primitive_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  Listing_primitive_traits(Output_iterator out_it, const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {}

  constexpr bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( m_traits.do_intersect_object()(query, primitive) )
    {
      *m_out_it++ = primitive.id();
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;
  const AABBTraits& m_traits;
};


/**
 * @class First_primitive_traits
 */
template<typename AABBTraits, typename Query>
class First_primitive_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  First_primitive_traits(const AABBTraits& traits)
    : m_is_found(false)
    , m_result()
    , m_traits(traits) {}

  bool go_further() const { return !m_is_found; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( m_traits.do_intersect_object()(query, primitive) )
    {
      m_result = boost::optional<typename Primitive::Id>(primitive.id());
      m_is_found = true;
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

  boost::optional<typename Primitive::Id> result() const { return m_result; }
  bool is_intersection_found() const { return m_is_found; }

private:
  bool m_is_found;
  boost::optional<typename Primitive::Id> m_result;
  const AABBTraits& m_traits;
};

/**
 * @class Do_intersect_traits
 */
template<typename AABBTraits, typename Query>
class Do_intersect_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  Do_intersect_traits(const AABBTraits& traits)
    : m_is_found(false), m_traits(traits)
  {}

  bool go_further() const { return !m_is_found; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( m_traits.do_intersect_object()(query, primitive) )
      m_is_found = true;
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

  bool is_intersection_found() const { return m_is_found; }

private:
  bool m_is_found;
  const AABBTraits& m_traits;
};


/**
 * @class Projection_traits
 */
template <typename AABBTraits>
class Projection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;

public:
  Projection_traits(const Point& hint,
                    const typename Primitive::Id& hint_primitive,
                    const AABBTraits& traits)
    : m_closest_point(hint),
      m_closest_primitive(hint_primitive),
      m_traits(traits)
  {}

  constexpr bool go_further() const { return true; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    Point new_closest_point = m_traits.closest_point_object()
      (query, primitive, m_closest_point);
    if( !m_traits.equal_3_object()(new_closest_point, m_closest_point) )
    {
      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point; // this effectively shrinks the sphere
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    return m_traits.compare_distance_object()
      (query, node.bbox(), m_closest_point) == CGAL::SMALLER;
  }

  Point closest_point() const { return m_closest_point; }
  Point_and_primitive_id closest_point_and_primitive() const
  {
    return Point_and_primitive_id(m_closest_point, m_closest_primitive);
  }

private:
  Point m_closest_point;
  typename Primitive::Id m_closest_primitive;
  const AABBTraits& m_traits;
};

}}} // end namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
