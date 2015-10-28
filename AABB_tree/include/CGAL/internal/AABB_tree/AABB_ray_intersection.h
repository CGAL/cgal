// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s) : Philipp Moeller
//

#ifndef CGAL_AABB_RAY_INTERSECTION_H
#define CGAL_AABB_RAY_INTERSECTION_H

#include <queue>
#include <functional>
#include <boost/optional.hpp>

namespace CGAL {

template<typename AABBTree>
class AABB_ray_intersection {
  typedef typename AABBTree::AABB_traits AABB_traits;
  typedef typename AABB_traits::Ray_3 Ray;
  typedef typename AABBTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
public:
  AABB_ray_intersection(const AABBTree& tree) : tree_(tree) {}

  boost::optional< Ray_intersection_and_primitive_id >
  ray_intersection(const Ray& query) const {
    // We hit the root, now continue on the children. Keep track of
    // nb_primitives through a variable in each Node on the stack. In
    // BVH_node::traversal this is done through the function parameter
    // nb_primitives in the recursion.
    typedef std::priority_queue< Node_ptr_with_ft, std::vector<Node_ptr_with_ft>,
                                 std::greater<Node_ptr_with_ft> > Heap_type;

    typename AABB_traits::Intersection
      intersection_obj = tree_.traits().intersection_object();
    typename AABB_traits::Intersection_distance
      intersection_distance_obj = tree_.traits().intersection_distance_object();
    Heap_type pq;
    boost::optional< Ray_intersection_and_primitive_id >
      intersection, /* the temporary for calculating the result */
      p; /* the current best intersection */

    // this is not the right way to do it, but using
    // numeric_limits<FT>::{max,infinity} will not work with Epeck.
    FT t = std::numeric_limits<double>::max();
    // Start with the root node.
    pq.push(Node_ptr_with_ft(tree_.root_node(), 0, tree_.size()));

    while(!pq.empty() && pq.top().value < t) {
      Node_ptr_with_ft current = pq.top();
      pq.pop();

      switch(current.nb_primitives) { // almost copy-paste from BVH_node::traversal
      case 2: // Left & right child both leaves
      {
        //left child
        intersection = intersection_obj(query, current.node->left_data());
        if(intersection) {
          FT ray_distance = as_ray_parameter(query, *intersection);
          if(ray_distance < t) {
            t = ray_distance;
            p = intersection;
          }
        }

        // right child
        intersection = intersection_obj(query, current.node->right_data());
        if(intersection) {
          FT ray_distance = as_ray_parameter(query, *intersection);
          if(ray_distance < t) {
            t = ray_distance;
            p = intersection;
          }
        }
      }
      case 3: // Left child leaf, right child inner node
      {
        //left child
        intersection = intersection_obj(query, current.node->left_data());
        if(intersection) {
          FT ray_distance = as_ray_parameter(query, *intersection);
          if(ray_distance < t) {
            t = ray_distance;
            p = intersection;
          }
        }

        // right child
        const Node* child = &(current.node->right_child());
        boost::optional< FT > dist = intersection_distance_obj(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, 2));

        break;
      }
      default: // Children both inner nodes
      {
        const Node* child = &(current.node->left_child());
        boost::optional<FT> dist = intersection_distance_obj(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, current.nb_primitives/2));

        child = &(current.node->right_child());
        dist = intersection_distance_obj(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, current.nb_primitives - current.nb_primitives/2));

        break;
      }
      }
    }

    return p;
  }
private:
  const AABBTree& tree_;
  typedef typename AABBTree::Point Point;
  typedef typename AABBTree::FT FT;
  typedef typename AABBTree::Node Node;
  typedef typename AABBTree::size_type size_type;

  struct Node_ptr_with_ft {
    Node_ptr_with_ft(const Node* node, const FT& value, size_type nb_primitives)
      : node(const_cast<Node*>(node)), value(value), nb_primitives(nb_primitives) {}
    Node* node;
    FT value;
    size_type nb_primitives;
    bool operator<(const Node_ptr_with_ft& other) const { return value < other.value; }
    bool operator>(const Node_ptr_with_ft& other) const { return value > other.value; }
  };

  FT as_ray_parameter(const Ray& ray,
                      const Ray_intersection_and_primitive_id& intersection)
    const {
    // TODO replace with non-hacky solution
    if(const Point* point = boost::get<const Point>(&(intersection.first))) {
      typename AABB_traits::Geom_traits::Vector_3 distance_ray(*point, ray.source());
      return distance_ray.squared_length();
    } else {
      std::cout << "not handled" << std::endl;
      return 0;
    }
  }
};

template<typename AABBTraits>
template<typename Ray>
boost::optional< typename AABB_tree<AABBTraits>::template Intersection_and_primitive_id<Ray>::Type >
AABB_tree<AABBTraits>::ray_intersection(const Ray& query) const {
  switch(size()) // copy-paste from AABB_tree::traversal
  {
  case 0: // Tree empty, nothing to intersect
    break;
  case 1: // Tree has 1 node, intersect directly
    return traits().intersection_object()(query, singleton_data());
  default: // Tree has >= 2 nodes
    if(traits().do_intersect_object()(query, root_node()->bbox())) {
      AABB_ray_intersection< AABB_tree<AABBTraits> > ri(*this);
      return ri.ray_intersection(query);
    } else {
      // but we don't hit the root
      break;
    }
  }
  return boost::none;
}

}


#endif /* CGAL_AABB_RAY_INTERSECTION_H */
