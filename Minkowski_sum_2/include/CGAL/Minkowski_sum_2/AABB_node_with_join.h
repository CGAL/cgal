// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_AABB_NODE_WITH_JOIN_H
#define CGAL_AABB_NODE_WITH_JOIN_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>

namespace CGAL {

/**
 * @class AABB_node_with_join
 *
 *
 */
template<typename AABBTraits>
class AABB_node_with_join
{
public:
  typedef typename AABBTraits::Bounding_box Bounding_box;

  /// Constructor
  AABB_node_with_join()
    : m_bbox()
    , m_p_left_child(NULL)
    , m_p_right_child(NULL)      { };

  /// Non virtual Destructor
  /// Do not delete children because the tree hosts and delete them
  ~AABB_node_with_join() { };

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

  /**
   * @param other_node root node of a tree which we want to traverse in parallel
   * @param traits the traversal traits that define the traversal behaviour
   * @param nb_primitives the number of primitives in this tree
   * @param nb_primitives_other the number of primitives in the other tree
   * @param first_stationary if true, the other_node is the translatable tree's root
   *
   * General traversal query for two trees.
   */
  template<class Traversal_traits>
  void traversal(const AABB_node_with_join &other_node,
                 Traversal_traits &traits,
                 const std::size_t nb_primitives,
                 const std::size_t nb_primitives_other,
                 bool first_stationary) const;

private:
  typedef AABBTraits AABB_traits;
  typedef AABB_node_with_join<AABB_traits> Node;
  typedef typename AABB_traits::Primitive Primitive;

  /// Helper functions
  const Node& left_child() const
                     { return *static_cast<Node*>(m_p_left_child); }
  const Node& right_child() const
                     { return *static_cast<Node*>(m_p_right_child); }
  const Primitive& left_data() const
                     { return *static_cast<Primitive*>(m_p_left_child); }
  const Primitive& right_data() const
                     { return *static_cast<Primitive*>(m_p_right_child); }

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
  typedef AABB_node_with_join<AABBTraits> Self;
  AABB_node_with_join(const Self& src);
  Self& operator=(const Self& src);

};  // end class AABB_node_with_join


template<typename Tr>
template<typename ConstPrimitiveIterator>
void
AABB_node_with_join<Tr>::expand(ConstPrimitiveIterator first,
                      ConstPrimitiveIterator beyond,
                      const std::size_t range,
                      const Tr& traits)
{
  m_bbox = traits.compute_bbox_object()(first, beyond);

  // sort primitives along longest axis aabb
  traits.sort_primitives_object()(first, beyond, m_bbox);

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
AABB_node_with_join<Tr>::traversal(const Query& query,
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

template<typename Tr>
template<class Traversal_traits>
void
AABB_node_with_join<Tr>::traversal(const AABB_node_with_join &other_node,
                         Traversal_traits &traits,
                         const std::size_t nb_primitives,
                         const std::size_t nb_primitives_other,
                         bool first_stationary) const
{
  if (nb_primitives >= nb_primitives_other)
  {
    switch(nb_primitives)
    {
      case 2: // Both trees contain 2 primitives, test all pairs
        traits.intersection(left_data(), other_node.left_data(), first_stationary);
        if (!traits.go_further()) return;
        traits.intersection(right_data(), other_node.right_data(), first_stationary);
        if (!traits.go_further()) return;
        traits.intersection(right_data(), other_node.left_data(), first_stationary);
        if (!traits.go_further()) return;
        traits.intersection(left_data(), other_node.right_data(), first_stationary);
        break;

      case 3: // This tree contains 3 primitives, the other 3 or 2
        // Both left children are primitives:
        traits.intersection(left_data(), other_node.left_data(), first_stationary);
        if (!traits.go_further()) return;

        // Test left child against all right leaves of the other tree
        if (nb_primitives_other == 2)
        {
          traits.intersection(left_data(), other_node.right_data(), first_stationary);
        }
        else
        {
          if (traits.do_intersect(left_data(), other_node.right_child(), first_stationary))
          {
            traits.intersection(left_data(), other_node.right_child().left_data(), first_stationary);
            if (!traits.go_further()) return;
            traits.intersection(left_data(), other_node.right_child().right_data(), first_stationary);
          }
        }
        if (!traits.go_further()) return;

        // Test right child against the other node
        if(traits.do_intersect(right_child(), other_node, first_stationary))
        {
          right_child().traversal(other_node, traits, 2, nb_primitives_other, first_stationary);
        }
        break;

      default: // This tree has two node-children, test both against the other node
        if( traits.do_intersect(left_child(), other_node, first_stationary) )
        {
          left_child().traversal(other_node, traits, nb_primitives/2, nb_primitives_other, first_stationary);
        }
        if (!traits.go_further()) return;
        if( traits.do_intersect(right_child(), other_node, first_stationary) )
        {
          right_child().traversal(other_node, traits, nb_primitives-nb_primitives/2, nb_primitives_other, first_stationary);
        }
    }
  }
  else
  {
    // The other node contains more primitives. Call this method the other way around:
    other_node.traversal(*this, traits, nb_primitives_other, nb_primitives, !first_stationary);
  }
}

} // end namespace CGAL

#endif // CGAL_AABB_NODE_WITH_JOIN_H
