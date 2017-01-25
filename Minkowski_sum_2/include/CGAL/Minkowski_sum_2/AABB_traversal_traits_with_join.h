// Copyright (c) 2008-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_TRAVERSAL_TRAITS_WITH_JOIN_H
#define CGAL_AABB_TRAVERSAL_TRAITS_WITH_JOIN_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/AABB_node_with_join.h>
#include <boost/optional.hpp>

namespace CGAL { 

namespace internal { namespace AABB_tree_with_join {

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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

public:
  typedef
  #if CGAL_INTERSECTION_VERSION < 2
  boost::optional<Object_and_primitive_id> 
  #else
  boost::optional< typename AABBTraits::template Intersection_and_primitive_id<Query>::Type >
  #endif
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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

public:
  Listing_intersection_traits(Output_iterator out_it, const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {}

  bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    #if CGAL_INTERSECTION_VERSION < 2
    boost::optional<Object_and_primitive_id>
    #else
    boost::optional< typename AABBTraits::template Intersection_and_primitive_id<Query>::Type >
    #endif
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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

public:
  Listing_primitive_traits(Output_iterator out_it, const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {}

  bool go_further() const { return true; }

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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

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
 * @class Do_intersect_joined_traits
 */
template<typename AABBTraits>
class Do_intersect_joined_traits
{
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef AABB_node_with_join<AABBTraits> Node;

public:

  Do_intersect_joined_traits(const Point &point) : m_is_found(false)
  {
    m_traits_ptr = new AABBTraits(point);
  }

  bool go_further() const { return !m_is_found; }

  void intersection(const Primitive &primitive1, const Primitive &primitive2, bool first_stationary)
  {
    if (first_stationary)
    {
      if (m_traits_ptr->do_intersect_object()(primitive1, primitive2))
      {
        m_is_found = true;
      }
    }
    else
    {
      if (m_traits_ptr->do_intersect_object()(primitive2, primitive1))
      {
        m_is_found = true;
      }
    }
  }

  bool do_intersect(const Node &node_1, const Node &node_2, bool first_stationary) const
  {
    if (first_stationary)
    {
      return m_traits_ptr->do_intersect_object()(node_1.bbox(), node_2.bbox());
    }
    else
    {
      return m_traits_ptr->do_intersect_object()(node_2.bbox(), node_1.bbox());
    }
  }

  bool do_intersect(const Node &node_1, const Primitive &primitive2, bool first_stationary) const
  {
    if (first_stationary)
    {
      return m_traits_ptr->do_intersect_object()(node_1.bbox(), primitive2);
    }
    else
    {
      return m_traits_ptr->do_intersect_object()(primitive2, node_1.bbox());
    }
  }

  bool do_intersect(const Primitive &primitive1, const Node &node_2, bool first_stationary) const
  {
    if (first_stationary)
    {
      return m_traits_ptr->do_intersect_object()(primitive1, node_2.bbox());
    }
    else
    {
      return m_traits_ptr->do_intersect_object()(node_2.bbox(), primitive1);
    }
  }

  bool is_intersection_found() const { return m_is_found; }

  ~Do_intersect_joined_traits() { delete m_traits_ptr; }

private:

  bool m_is_found;
  AABBTraits *m_traits_ptr;
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
  typedef ::CGAL::AABB_node_with_join<AABBTraits> Node;

public:
  Projection_traits(const Point& hint,
                    const typename Primitive::Id& hint_primitive,
                    const AABBTraits& traits)
    : m_closest_point(hint),
      m_closest_primitive(hint_primitive), 
      m_traits(traits)
  {}

  bool go_further() const { return true; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    Point new_closest_point = m_traits.closest_point_object()
      (query, primitive, m_closest_point);
    if(new_closest_point != m_closest_point)
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

}}} // end namespace CGAL::internal::AABB_tree_with_join

#endif // CGAL_AABB_TRAVERSAL_TRAITS_WITH_JOIN_H
