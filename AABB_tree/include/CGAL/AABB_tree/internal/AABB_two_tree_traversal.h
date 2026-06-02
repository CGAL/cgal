// Copyright (c) 2026  Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_AABB_TWO_TREE_TRAVERSAL_TRAITS_H
#define CGAL_AABB_TWO_TREE_TRAVERSAL_TRAITS_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <optional>

namespace CGAL {

namespace internal { namespace AABB_tree {

template< typename Value, typename OutputIterator >
class WrapOutputIterator
{
  Value first;
  OutputIterator out;
  bool in_order;
public:
  WrapOutputIterator(Value first_, OutputIterator out_, bool in_order_ = true): first(first_), out(out_), in_order(in_order_){}

  WrapOutputIterator& operator=(Value second){
    if (in_order)
      *out = std::make_pair(first, second);
    else
      *out = std::make_pair(second, first);
    return *this;
  }

  WrapOutputIterator& operator*(){ return *this; }
  WrapOutputIterator& operator++(){ ++out; return *this; }
  WrapOutputIterator& operator++(int){ ++out; return *this; }
  WrapOutputIterator& operator+(int d){ auto tmp = *this; ++out; return *this; }
};

template <bool in_order = true,
          class Two_tree_traversal_traits,
          class AABBTraits_A,
          class AABBTraits_B,
          class OutputIterator,
          class Bbox_map_A,
          class Bbox_map_B>
void two_tree_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                        const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                        const std::size_t nb_primitives_A,
                        const std::size_t nb_primitives_B,
                        const AABBTraits_A& traits_A,
                        const AABBTraits_B& traits_B,
                        const Bbox_map_A &bb_A,
                        const Bbox_map_B &bb_B,
                        OutputIterator out)
{
  using Primitive_id = typename AABBTraits_A::Primitive::Id;
  using WrapOutputIterator = WrapOutputIterator<Primitive_id, OutputIterator>;
  using Intersection_traits = CGAL::internal::AABB_tree::Listing_primitive_traits<AABBTraits_B, Bbox_3, WrapOutputIterator>;

  switch(nb_primitives_A)
  {
  case 2:
  {
    WrapOutputIterator wrap_out_left(node_A.left_data().id(), out, in_order);
    Two_tree_traversal_traits traits_left(wrap_out_left, traits_B);
    node_B.traversal( get(bb_A, node_A.left_data().id()), traits_left, nb_primitives_B);

    WrapOutputIterator wrap_out_right(node_A.right_data().id(), out, in_order);
    Two_tree_traversal_traits traits_right(wrap_out_right, traits_B);
    node_B.traversal( get(bb_A, node_A.right_data().id()) , traits_right, nb_primitives_B);
    break;
  }
  case 3:
  {
    WrapOutputIterator wrap_out_left(node_A.left_data().id(), out, in_order);
    Two_tree_traversal_traits traits_left(wrap_out_left, traits_B);
    node_B.traversal( get(bb_A, node_A.left_data().id()), traits_left, nb_primitives_B);

    if( do_overlap(node_A.right_child().bbox(), node_B.bbox()) )
      two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.right_child(), nb_primitives_B, 2, traits_B, traits_A, bb_B, bb_A, out);
    break;
  }
  default:
  {
    // TODO if both intersect, use parallel_invoke
    // TODO current strategy is we swap at each step. Try largest tree first or largest bbox first
    bool do_intersect_left  = do_overlap(node_A.left_child().bbox(), node_B.bbox());
    bool do_intersect_right = do_overlap(node_A.right_child().bbox(), node_B.bbox());
    if( do_intersect_left )
      two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traits_B, traits_A, bb_B, bb_A, out);
    if( do_intersect_right )
      two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2,  traits_B, traits_A, bb_B, bb_A, out);
  }} // switch end
}

template<class Tree_A,
         class Tree_B,
         class OutputIterator,
         class Bbox_map_A,
         class Bbox_map_B>
void two_tree_listing_intersecting_primitives(const Tree_A& tree_A,
                                              const Tree_B& tree_B,
                                              OutputIterator out,
                                              const Bbox_map_A &bb_A,
                                              const Bbox_map_B &bb_B)
{
  using Primitive_id = typename Tree_A::AABB_traits::Primitive::Id;
  using WrapOutputIterator = WrapOutputIterator<Primitive_id, OutputIterator>;
  using Two_tree_Listing_primitive_traits = CGAL::internal::AABB_tree::Listing_primitive_traits<typename Tree_B::AABB_traits, Bbox_3, WrapOutputIterator>;
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_tree_traversal<true, Two_tree_Listing_primitive_traits>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), tree_A.traits(), tree_B.traits(), bb_A, bb_B, out);
}

}}} // end namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
