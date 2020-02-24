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
public:
  typedef typename AABBTraits::Bounding_box Bounding_box;

  /// Constructor
  AABB_node()
    : m_bbox()
    , m_p_left_child(nullptr)
    , m_p_right_child(nullptr)      { };

  /// Non virtual Destructor
  /// Do not delete children because the tree hosts and delete them
  ~AABB_node() { };

  /// Returns the bounding box of the node
  const Bounding_box& bbox() const { return m_bbox; }

  /**
   * @brief Builds the tree by recursive expansion.
   * @param first the first primitive to insert
   * @param last the last primitive to insert
   * @param range the number of primitive of the range
   *
   * [first,last[ is the range of primitives to be added to the tree.
   */
  template<typename ConstPrimitiveIterator>
  void expand(ConstPrimitiveIterator first,
              ConstPrimitiveIterator beyond,
              const std::size_t range,
              const AABBTraits&);

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

private:
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

private:
  // Disabled copy constructor & assignment operator
  typedef AABB_node<AABBTraits> Self;
  AABB_node(const Self& src);
  Self& operator=(const Self& src);

};  // end class AABB_node


template<typename Tr>
template<typename ConstPrimitiveIterator>
void
AABB_node<Tr>::expand(ConstPrimitiveIterator first,
                      ConstPrimitiveIterator beyond,
                      const std::size_t range,
                      const Tr& traits)
{
  m_bbox = traits.compute_bbox_object()(first, beyond);

  // sort primitives along longest axis aabb
  traits.split_primitives_object()(first, beyond, m_bbox);

  switch(range)
  {
  case 2:
    m_p_left_child = &(*first);
    m_p_right_child = &(*(++first));
    break;
  case 3:
    m_p_left_child = &(*first);
    m_p_right_child = static_cast<Node*>(this)+1;
    right_child().expand(first+1, beyond, 2,traits);
    break;
  default:
    const std::size_t new_range = range/2;
    m_p_left_child = static_cast<Node*>(this) + 1;
    m_p_right_child = static_cast<Node*>(this) + new_range;
    left_child().expand(first, first + new_range, new_range,traits);
    right_child().expand(first + new_range, beyond, range - new_range,traits);
  }
}


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
