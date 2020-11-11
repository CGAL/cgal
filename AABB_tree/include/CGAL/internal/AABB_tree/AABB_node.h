// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_AABB_NODE_H
#define CGAL_AABB_NODE_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>

namespace CGAL {

/**
 * @class AABB_node
 *
 *
 */
template<typename AABBTraits>
class AABB_node
{
private:
  typedef AABB_node<AABBTraits> Self;

public:
  typedef typename AABBTraits::Bounding_box Bounding_box;

  /// Constructor
  AABB_node()
    : m_bbox()
    , m_p_left_child(nullptr)
    , m_p_right_child(nullptr)      { };

  AABB_node(Self&& node) = default;

  // Disabled copy constructor & assignment operator
  AABB_node(const Self& src) = delete;
  Self& operator=(const Self& src) = delete;

  /// Returns the bounding box of the node
  const Bounding_box& bbox() const { return m_bbox; }

  /**
   * @brief General traversal query
   * @param query the query
   * @param traits the traversal traits that define the traversal behaviour
   * @param nb_primitives the number of primitive
   *
   * General traversal query. The traits class allows using it for the various
   * traversal methods we need: listing, counting, detecting intersections,
   * drawing the boxes.
   */
  template<class Traversal_traits, class Query>
  void traversal(const Query& query,
                 Traversal_traits& traits,
                 const std::size_t nb_primitives) const;

private:
  typedef AABBTraits AABB_traits;
  typedef AABB_node<AABB_traits> Node;
  typedef typename AABB_traits::Primitive Primitive;


public:
  /// Helper functions
  const Node& left_child() const
                     { return *static_cast<Node*>(m_p_left_child); }
  const Node& right_child() const
                     { return *static_cast<Node*>(m_p_right_child); }
  const Primitive& left_data() const
                     { return *static_cast<Primitive*>(m_p_left_child); }
  const Primitive& right_data() const
                     { return *static_cast<Primitive*>(m_p_right_child); }
  template <class Left, class Right>
  void set_children(Left& l, Right& r)
  {
    m_p_left_child = static_cast<void*>(std::addressof(l));
    m_p_right_child = static_cast<void*>(std::addressof(r));
  }
  void set_bbox(const Bounding_box& bbox)
  {
    m_bbox = bbox;
  }

  Node& left_child() { return *static_cast<Node*>(m_p_left_child); }
  Node& right_child() { return *static_cast<Node*>(m_p_right_child); }
  Primitive& left_data() { return *static_cast<Primitive*>(m_p_left_child); }
  Primitive& right_data() { return *static_cast<Primitive*>(m_p_right_child); }

private:
  /// node bounding box
  Bounding_box m_bbox;

  /// children nodes, either pointing towards children (if children are not leaves),
  /// or pointing toward input primitives (if children are leaves).
  void *m_p_left_child;
  void *m_p_right_child;

};  // end class AABB_node


template<typename Tr>
template<class Traversal_traits, class Query>
void
AABB_node<Tr>::traversal(const Query& query,
                         Traversal_traits& traits,
                         const std::size_t nb_primitives) const
{
  // Recursive traversal
  switch(nb_primitives)
  {
  case 2:
    traits.intersection(query, left_data());
    if( traits.go_further() )
    {
      traits.intersection(query, right_data());
    }
    break;
  case 3:
    traits.intersection(query, left_data());
    if( traits.go_further() && traits.do_intersect(query, right_child()) )
    {
      right_child().traversal(query, traits, 2);
    }
    break;
  default:
    if( traits.do_intersect(query, left_child()) )
    {
      left_child().traversal(query, traits, nb_primitives/2);
      if( traits.go_further() && traits.do_intersect(query, right_child()) )
      {
        right_child().traversal(query, traits, nb_primitives-nb_primitives/2);
      }
    }
    else if( traits.do_intersect(query, right_child()) )
    {
      right_child().traversal(query, traits, nb_primitives-nb_primitives/2);
    }
  }
}

} // end namespace CGAL

#endif // CGAL_AABB_NODE_H
