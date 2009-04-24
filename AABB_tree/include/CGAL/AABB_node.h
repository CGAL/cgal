// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// Copyrigth (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     :  Camille Wormser, Jane Tournois, Pierre Alliez, Laurent Rineau

#ifndef CGAL_AABB_NODE_H
#define CGAL_AABB_NODE_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>

#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>

#include <CGAL/AABB_intersections.h>

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
    : bbox_()
    , p_left_child_(NULL)
    , p_right_child_(NULL)      { };

  /// Non virtual Destructor
  /// Do not delete children because the tree hosts and delete them
  ~AABB_node() { };

  /// Returns the bounding box of the node
  Bounding_box bounding_box() const { return bbox_; }

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
              ConstPrimitiveIterator last,
              const int range);


  /**
   * @brief General traversal query
   * @param query the query
   * @param traits the traversal traits that define the traversal behaviour
   * @param nb_primitives the number of primitive
   *
   * General traversal query, the traits class allows to use it for the various
   * traversal methods we need: listing, counting, detecting intersections,
   * drawing the boxes...
   */
  template<class Traversal_traits, class Query>
  void traversal(const Query& query,
                 Traversal_traits& traits,
                 const int nb_primitives) const;


private:
  typedef AABB_node<AABBTraits> Node;
  typedef typename AABBTraits::Primitive Primitive;



  /// Helper functions
  const Node& left_child() const
                     { return *static_cast<Node*>(p_left_child_); }
  const Node& right_child() const
                     { return *static_cast<Node*>(p_right_child_); }
  const Primitive& left_data() const
                     { return *static_cast<Primitive*>(p_left_child_); }
  const Primitive& right_data() const
                     { return *static_cast<Primitive*>(p_right_child_); }

  Node& left_child() { return *static_cast<Node*>(p_left_child_); }
  Node& right_child() { return *static_cast<Node*>(p_right_child_); }
  Primitive& left_data() { return *static_cast<Primitive*>(p_left_child_); }
  Primitive& right_data() { return *static_cast<Primitive*>(p_right_child_); }

private:
  /// bounding box
  Bounding_box bbox_;

  /// children nodes:
  /// either pointing towards children (if the children are not leaves)
  /// or pointing toward input primitives (if the children are leaves).
  void *p_left_child_;
  void *p_right_child_;

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
                      ConstPrimitiveIterator last,
                      const int range)
{
  bbox_ = Tr().compute_bbox(first, last);

  // sort primitives along longest axis aabb
  Tr().sort_primitives(first, last, bbox_);

  switch(range)
  {
  case 2:
    p_left_child_ = &(*first);
    p_right_child_ = &(*(++first));
    break;
  case 3:
    p_left_child_ = &(*first);
    p_right_child_ = static_cast<Node*>(this)+1;
    right_child().expand(first+1, last, 2);
    break;
  default:
    const int new_range = range/2;
    p_left_child_ = static_cast<Node*>(this)+1;
    p_right_child_ = static_cast<Node*>(this)+new_range;
    left_child().expand(first, first+new_range, new_range);
    right_child().expand(first+new_range, last, range-new_range);
  }
}


template<typename Tr>
template<class Traversal_traits, class Query>
void
AABB_node<Tr>::traversal(const Query& query,
                         Traversal_traits& traits,
                         const int nb_primitives) const
{
//  CGAL_assertion(NULL!=p_left_child_);
//  CGAL_assertion(NULL!=p_right_child_);

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

