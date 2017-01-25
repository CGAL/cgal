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

#include <CGAL/license/AABB_tree.h>


#include <functional>
#include <boost/optional.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/variant/apply_visitor.hpp>
#if BOOST_VERSION >= 105000
#include <boost/heap/priority_queue.hpp>
#else
#include <queue>
#endif
#include <CGAL/assertions.h>

namespace CGAL {

template<typename AABBTree, typename SkipFunctor>
class AABB_ray_intersection {
  typedef typename AABBTree::AABB_traits AABB_traits;
  typedef typename AABB_traits::Ray_3 Ray;
  typedef typename AABBTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
  typedef typename Ray_intersection_and_primitive_id::first_type Ray_intersection;
public:
  AABB_ray_intersection(const AABBTree& tree) : tree_(tree) {}

  boost::optional< Ray_intersection_and_primitive_id >
  ray_intersection(const Ray& query, SkipFunctor skip) const {
    // We hit the root, now continue on the children. Keep track of
    // nb_primitives through a variable in each Node on the stack. In
    // BVH_node::traversal this is done through the function parameter
    // nb_primitives in the recursion.
        typedef
#if BOOST_VERSION >= 105000
          boost::heap::priority_queue< Node_ptr_with_ft, boost::heap::compare< std::greater<Node_ptr_with_ft> > >
#else
          std::priority_queue< Node_ptr_with_ft>
#endif
          Heap_type;

    typename AABB_traits::Intersection
      intersection_obj = tree_.traits().intersection_object();
    // typename AABB_traits::Do_intersect
    //   do_intersect_obj = tree_.traits().do_intersect_object();
    typename AABB_traits::Intersection_distance
      intersection_distance_obj = tree_.traits().intersection_distance_object();
    as_ray_param_visitor param_visitor = as_ray_param_visitor(&query);

    Heap_type pq;
    // pq.reserve(tree_.size() / 2);
    boost::optional< Ray_intersection_and_primitive_id >
      intersection, /* the temporary for calculating the result */
      p; /* the current best intersection */

    // this is not the right way to do it, but using
    // numeric_limits<FT>::{max,infinity} will not work with Epeck.
    FT t = (std::numeric_limits<double>::max)();
    // Start with the root node.
    pq.push(Node_ptr_with_ft(tree_.root_node(), 0, tree_.size()));

    while(!pq.empty() && pq.top().value < t) {
      Node_ptr_with_ft current = pq.top();
      pq.pop();

      switch(current.nb_primitives) { // almost copy-paste from BVH_node::traversal
      case 2: // Left & right child both leaves
      {
        //left child
        if(!skip(current.node->left_data().id()) /* && do_intersect_obj(query, current.node->left_data()) */) {
          intersection = intersection_obj(query, current.node->left_data());
          if(intersection) {
            FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
            if(ray_distance < t) {
              t = ray_distance;
              p = intersection;
            }
          }
        }

        // right child
        if(!skip(current.node->right_data().id()) /* && do_intersect_obj(query, current.node->right_data()) */) {
          intersection = intersection_obj(query, current.node->right_data());
          if(intersection) {
            FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
            if(ray_distance < t) {
              t = ray_distance;
              p = intersection;
            }
          }
        }
        break;
      }
      case 3: // Left child leaf, right child inner node
      {
        //left child
        if(!skip(current.node->left_data().id()) /* && do_intersect_obj(query, current.node->left_data()) */) {
          intersection = intersection_obj(query, current.node->left_data());
          if(intersection) {
            FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
            if(ray_distance < t) {
              t = ray_distance;
              p = intersection;
            }
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
      : node(node), nb_primitives(nb_primitives), value(value) {}
    const Node* node;
    size_type nb_primitives;
    FT value;
#if BOOST_VERSION >= 105000
    bool operator<(const Node_ptr_with_ft& other) const { return value < other.value; }
    bool operator>(const Node_ptr_with_ft& other) const { return value > other.value; }
#else
    bool operator>(const Node_ptr_with_ft& other) const { return value < other.value; }
    bool operator<(const Node_ptr_with_ft& other) const { return value > other.value; }

#endif
  };

  struct as_ray_param_visitor {
    typedef FT result_type;
    as_ray_param_visitor(const Ray* ray) : ray(ray) {}

    template<typename T>
    FT operator()(const T& s)
    {
      // intersection is a segment, returns the min relative distance
      // of its endpoints
      FT r1 = this->operator()(s[0]);
      FT r2 = this->operator()(s[1]);
      return (std::min)(r1,r2);
    }

    FT operator()(const Point& point) {
      typename AABB_traits::Geom_traits::Vector_3 x(ray->source(), point);
      typename AABB_traits::Geom_traits::Vector_3 v = ray->to_vector();

      for(int i = 0; i < 3; ++i) {
        if(v[i] != FT(0.)) {
          return x[i] / v[i];
        }
      }
      CGAL_assertion(false); // should never end-up here
      return FT(0.);
    }

    const Ray* ray;
  };
};

template<typename AABBTraits>
template<typename Ray, typename SkipFunctor>
boost::optional< typename AABB_tree<AABBTraits>::template Intersection_and_primitive_id<Ray>::Type >
AABB_tree<AABBTraits>::first_intersection(const Ray& query,
                                          const SkipFunctor& skip) const {
  CGAL_static_assertion_msg((boost::is_same<Ray, typename AABBTraits::Ray_3>::value), 
                            "Ray and Ray_3 must be the same type");

  switch(size()) // copy-paste from AABB_tree::traversal
  {
  case 0: // Tree empty, nothing to intersect
    break;
  case 1: // Tree has 1 node, intersect directly
    return traits().intersection_object()(query, singleton_data());
  default: // Tree has >= 2 nodes
    if(traits().do_intersect_object()(query, root_node()->bbox())) {
      AABB_ray_intersection< AABB_tree<AABBTraits>, SkipFunctor > ri(*this);
      return ri.ray_intersection(query, skip);
    } else {
      // but we don't hit the root
      break;
    }
  }
  return boost::none;
}

template<typename AABBTraits>
template<typename Ray, typename SkipFunctor>
boost::optional<typename AABB_tree<AABBTraits>::Primitive_id>
AABB_tree<AABBTraits>::first_intersected_primitive(const Ray& query,
                            const SkipFunctor& skip) const
{
  boost::optional<
    typename AABB_tree<AABBTraits>::
      template Intersection_and_primitive_id<Ray>::Type > res =
        first_intersection(query, skip);
  if ( (bool) res )
    return boost::make_optional( res->second );
  return boost::none;
}

}

#endif /* CGAL_AABB_RAY_INTERSECTION_H */
