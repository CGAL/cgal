// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Philipp Moeller
//

#ifndef CGAL_AABB_RAY_INTERSECTION_H
#define CGAL_AABB_RAY_INTERSECTION_H

#include <CGAL/license/AABB_tree.h>


#include <functional>
#include <type_traits>
#include <optional>
#  if defined(BOOST_MSVC)
#    pragma warning(push)
#    pragma warning(disable: 4996)
#  endif
#  include <boost/heap/priority_queue.hpp>
#  if defined(BOOST_MSVC)
#    pragma warning(pop)
#  endif

#include <CGAL/assertions.h>

namespace CGAL {

template<typename AABBTree, typename SkipFunctor>
class AABB_ray_intersection {
  typedef typename AABBTree::AABB_traits AABB_traits;
  static const int dimension = AABB_traits::Point::Ambient_dimension::value;
  typedef typename AABB_traits::Ray Ray;
  typedef typename AABB_traits::Vector Vector;

  typedef typename AABBTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
  typedef typename Ray_intersection_and_primitive_id::first_type Ray_intersection;

public:
  AABB_ray_intersection(const AABBTree& tree) : tree_(tree) {}

  std::optional< Ray_intersection_and_primitive_id >
  ray_intersection(const Ray& query, SkipFunctor skip) const {
    // We hit the root, now continue on the children. Keep track of
    // nb_primitives through a variable in each Node on the stack. In
    // BVH_node::traversal this is done through the function parameter
    // nb_primitives in the recursion.
    typedef boost::heap::priority_queue< Node_ptr_with_ft, boost::heap::compare< std::greater<Node_ptr_with_ft> > >
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
    std::optional< Ray_intersection_and_primitive_id >
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
            FT ray_distance = std::visit(param_visitor, intersection->first);
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
            FT ray_distance = std::visit(param_visitor, intersection->first);
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
            FT ray_distance = std::visit(param_visitor, intersection->first);
            if(ray_distance < t) {
              t = ray_distance;
              p = intersection;
            }
          }
        }

        // right child
        const Node* child = &(current.node->right_child());
        std::optional< FT > dist = intersection_distance_obj(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, 2));

        break;
      }
      default: // Children both inner nodes
      {
        const Node* child = &(current.node->left_child());
        std::optional<FT> dist = intersection_distance_obj(query, child->bbox());
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
    bool operator<(const Node_ptr_with_ft& other) const { return value < other.value; }
    bool operator>(const Node_ptr_with_ft& other) const { return value > other.value; }
  };

  struct as_ray_param_visitor {
    typedef FT result_type;
    as_ray_param_visitor(const Ray* ray)
     : ray(ray), max_i(0)
    {
      Vector v = AABB_traits().construct_vector_object()(*ray);
      for (int i=1; i<dimension; ++i)
        if( CGAL::abs(v[i]) > CGAL::abs(v[max_i]) )
          max_i = i;
    }

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
      Vector x = Vector(AABB_traits().construct_source_object()(*ray), point);
      Vector v = AABB_traits().construct_vector_object()(*ray);

      return x[max_i] / v[max_i];
    }

    const Ray* ray;
    int max_i;
  };
};

template<typename AABBTraits>
template<typename Ray, typename SkipFunctor>
std::optional< typename AABB_tree<AABBTraits>::template Intersection_and_primitive_id<Ray>::Type >
AABB_tree<AABBTraits>::first_intersection(const Ray& query,
                                          const SkipFunctor& skip) const {
  static_assert(std::is_same<Ray, typename AABBTraits::Ray>::value,
                            "Ray and AABBTraits::Ray must be the same type");

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
  return std::nullopt;
}

template<typename AABBTraits>
template<typename Ray, typename SkipFunctor>
std::optional<typename AABB_tree<AABBTraits>::Primitive_id>
AABB_tree<AABBTraits>::first_intersected_primitive(const Ray& query,
                            const SkipFunctor& skip) const
{
  std::optional<
    typename AABB_tree<AABBTraits>::
      template Intersection_and_primitive_id<Ray>::Type > res =
        first_intersection(query, skip);
  if ( (bool) res )
    return std::make_optional( res->second );
  return std::nullopt;
}

}

#endif /* CGAL_AABB_RAY_INTERSECTION_H */
