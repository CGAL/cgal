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

#ifndef CGAL_AABB_TRAVERSAL_TRAITS_H
#define CGAL_AABB_TRAVERSAL_TRAITS_H

#include <CGAL/internal/AABB_tree/AABB_node.h>
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
  typedef typename boost::optional<Object_and_primitive_id> Result;
public:
  First_intersection_traits()
    : m_result()
  {}

  bool go_further() const { return !m_result; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    m_result = AABBTraits().intersection_object()(query, primitive);
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

  Result result() const { return m_result; }
  bool is_intersection_found() const { return m_result; }

private:
  Result m_result;
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
  Listing_intersection_traits(Output_iterator out_it)
    : m_out_it(out_it) {}

  bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    boost::optional<Object_and_primitive_id> intersection;
    intersection = AABBTraits().intersection_object()(query, primitive);
    if(intersection)
    {
      *m_out_it++ = *intersection;
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;
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
  Listing_primitive_traits(Output_iterator out_it)
    : m_out_it(out_it) {}

  bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( AABBTraits().do_intersect_object()(query, primitive) )
    {
      *m_out_it++ = primitive.id();
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;
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
  First_primitive_traits()
    : m_is_found(false)
    , m_result() {}

  bool go_further() const { return !m_is_found; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( AABBTraits().do_intersect_object()(query, primitive) )
    {
      m_result = boost::optional<typename Primitive::Id>(primitive.id());
      m_is_found = true;
    }
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

  boost::optional<typename Primitive::Id> result() const { return m_result; }
  bool is_intersection_found() const { return m_is_found; }

private:
  bool m_is_found;
  boost::optional<typename Primitive::Id> m_result;
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
  Do_intersect_traits()
    : m_is_found(false)
  {}

  bool go_further() const { return !m_is_found; }

  void intersection(const Query& query, const Primitive& primitive)
  {
    if( AABBTraits().do_intersect_object()(query, primitive) )
      m_is_found = true;
  }

  bool do_intersect(const Query& query, const Node& node) const
  {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

  bool is_intersection_found() const { return m_is_found; }

private:
  bool m_is_found;
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
                    const typename Primitive::Id& hint_primitive)
    : m_closest_point(hint),
      m_closest_primitive(hint_primitive)
  {}

  bool go_further() const { return true; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    Point new_closest_point = AABBTraits().closest_point_object()
      (query, primitive, m_closest_point);
    if(new_closest_point != m_closest_point)
    {
      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point; // this effectively shrinks the sphere 
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    return AABBTraits().compare_distance_object()
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
};

}}} // end namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
