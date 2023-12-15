// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/offset_intersection.h>

#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <CGAL/Default.h>

#include <algorithm>
#include <memory>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename AABBTraits>
struct Default_traversal_traits
{
  using Projection_traits = CGAL::internal::AABB_tree::Projection_traits<AABBTraits>;

  template <typename Query>
  using Do_intersect_traits = CGAL::internal::AABB_tree::Do_intersect_traits<AABBTraits, Query>;

  template <typename Query>
  using First_intersection_traits = CGAL::internal::AABB_tree::First_intersection_traits<AABBTraits, Query>;
};

// Factorize the implementation of the functions calling the AABB tree
template <typename AABBTree,
          typename AABBTraversalTraits>
struct AABB_tree_oracle_helper
{
  using Self = AABB_tree_oracle_helper<AABBTree, AABBTraversalTraits>;

  using AABB_traits = typename AABBTree::AABB_traits;
  using GT = typename AABB_traits::Geom_traits;

  using FT = typename AABB_traits::FT;
  using Point_3 = typename AABB_traits::Point_3;

  template <typename Query>
  static bool do_intersect(const Query& query,
                           const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using Do_intersect_traits = typename AABBTraversalTraits::template Do_intersect_traits<Query>;

    Do_intersect_traits traversal_traits(tree.traits());
    tree.traversal(query, traversal_traits);
    return traversal_traits.is_intersection_found();
  }

  static Point_3 closest_point(const Point_3& p,
                               const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using Projection_traits = typename AABBTraversalTraits::Projection_traits;

    const auto& hint = tree.best_hint(p);

    Projection_traits projection_traits(hint.first, hint.second, tree.traits());
    tree.traversal(p, projection_traits);
    return projection_traits.closest_point();
  }

  static FT squared_distance(const Point_3& p,
                             const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    const Point_3 closest = Self::closest_point(p, tree);
    return tree.traits().squared_distance_object()(p, closest);
  }

  static bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o,
                                 const FT offset_size,
                                 const FT intersection_precision,
                                 const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using AABB_distance_oracle = internal::AABB_distance_oracle<AABBTree, AABBTraversalTraits>;
    using Offset_intersection = internal::Offset_intersection<GT, AABB_distance_oracle>;

    AABB_distance_oracle dist_oracle(tree);
    Offset_intersection offset_intersection(dist_oracle, offset_size, intersection_precision, 1 /*lip*/);
    return offset_intersection.first_intersection(p, q, o);
  }
};

template <typename GT,
          typename AABBTree,
          typename AABBTraversalTraits = CGAL::Default,
          typename BaseOracle = int> // base oracle
class AABB_tree_oracle
  : public BaseOracle
{
protected:
  using Geom_traits = GT;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  // When building oracle stacks, there are copies of (empty) trees, which isn't possible, thus pointer
  using AABB_tree = AABBTree;
  using AABB_tree_ptr = std::shared_ptr<AABB_tree>;
  using AABB_traits = typename AABB_tree::AABB_traits;
  using AABB_traversal_traits = typename Default::Get<AABBTraversalTraits,
                                                      Default_traversal_traits<AABB_traits> >::type;

  using AABB_helper = AABB_tree_oracle_helper<AABB_tree, AABB_traversal_traits>;

public:
  AABB_tree_oracle(const BaseOracle& base,
                   const GT& gt)
    : BaseOracle(base),
      m_gt(gt),
      m_tree_ptr(std::make_shared<AABB_tree>())
  { }

  AABB_tree_oracle(const AABB_tree_oracle&) = default;

public:
  const Geom_traits& geom_traits() const { return m_gt; }

  AABB_tree& tree() { return *m_tree_ptr; }
  const AABB_tree& tree() const { return *m_tree_ptr; }
  BaseOracle& base() { return static_cast<BaseOracle&>(*this); }
  const BaseOracle& base() const { return static_cast<const BaseOracle&>(*this); }

  bool empty() const { return m_tree_ptr->empty(); }
  bool do_call() const { return (!empty() || base().do_call()); }

  void clear() { m_tree_ptr->clear() && base().clear(); }

public:
  typename AABB_tree::Bounding_box bbox() const
  {
    CGAL_precondition(do_call());

    typename AABB_tree::Bounding_box bb;

    if(base().do_call())
      bb += base().bbox();

    if(!empty())
      bb += tree().bbox();

    return bb;
  }

  template <typename T>
  bool do_intersect(const T& t) const
  {
    CGAL_precondition(do_call());

    if(base().do_call() && base().do_intersect(t))
      return true;

    if(!empty())
      return AABB_helper::do_intersect(t, tree());

    return false;
  }

  FT squared_distance(const Point_3& p) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        const FT base_sqd = base().squared_distance(p);
        // @speed could do a smarter traversal, no need to search deeper than the current best
        const FT this_sqd = AABB_helper::squared_distance(p, tree());
        return (std::min)(base_sqd, this_sqd);
      }
      else // this level is empty
      {
        return base().squared_distance(p);
      }
    }
    else // empty base
    {
      return AABB_helper::squared_distance(p, tree());
    }
  }

  Point_3 closest_point(const Point_3& p) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        const Point_3 base_c = base().closest_point(p);
        // @speed could do a smarter traversal, no need to search deeper than the current best
        const Point_3 this_c = AABB_helper::closest_point(p, tree());
        return (compare_distance_to_point(p, base_c, this_c) == CGAL::SMALLER) ? base_c : this_c;
      }
      else // this level is empty
      {
        return base().closest_point(p);
      }
    }
    else // empty base
    {
      return AABB_helper::closest_point(p, tree());
    }
  }

  bool first_intersection(const Point_3& p, const Point_3& q,
                          Point_3& o,
                          const FT offset_size,
                          const FT intersection_precision) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        Point_3 base_o;
        bool base_b = base().first_intersection(p, q, base_o, offset_size, intersection_precision);

        if(base_b) // intersection found in base
        {
          // @speed could do a smarter traversal, no need to search deeper than the current best
          Point_3 this_o;
          bool this_b = AABB_helper::first_intersection(p, q, this_o, offset_size, intersection_precision, tree());
          if(this_b)
            o = (compare_distance_to_point(p, base_o, this_o) == SMALLER) ? base_o : this_o;
          else
            o = base_o;

          return true;
        }
        else // no intersection found in non-empty base
        {
          return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
        }
      }
      else // this level is empty
      {
        return base().first_intersection(p, q, o, offset_size, intersection_precision);
      }
    }
    else // empty base
    {
      return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
    }
  }

  bool first_intersection(const Point_3& p, const Point_3& q,
                          Point_3& o,
                          const FT offset_size) const
  {
    return first_intersection(p, q, o, offset_size, 1e-2 * offset_size);
  }

protected:
  Geom_traits m_gt;
  AABB_tree_ptr m_tree_ptr;
};

// partial specialization, when there is no further oracle beneath in the stack.
//
// `int` is used to denote the absence of base rather than `void`,
// as to use the same constructor for both versions (thus requires a default construction)
template <typename GT,
          typename AABBTree,
          typename AABBTraversalTraits>
class AABB_tree_oracle<GT, AABBTree, AABBTraversalTraits, int>
{
protected:
  using Geom_traits = GT;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  using AABB_tree = AABBTree;
  using AABB_tree_ptr = std::shared_ptr<AABB_tree>;
  using AABB_traits = typename AABB_tree::AABB_traits;
  using AABB_traversal_traits = typename Default::Get<AABBTraversalTraits,
                                                      Default_traversal_traits<AABB_traits> >::type;

  using AABB_helper = AABB_tree_oracle_helper<AABB_tree, AABB_traversal_traits>;

public:
  AABB_tree_oracle(const int, // to have a common constructor API between both versions
                   const GT& gt)
    : m_gt(gt), m_tree_ptr(std::make_shared<AABB_tree>())
  { }

public:
  const Geom_traits& geom_traits() const { return m_gt; }
  AABB_tree& tree() { return *m_tree_ptr; }
  const AABB_tree& tree() const { return *m_tree_ptr; }

  bool empty() const { return m_tree_ptr->empty(); }
  bool do_call() const { return !empty(); }

  void clear() { m_tree_ptr->clear(); }

public:
  typename AABB_tree::Bounding_box bbox() const
  {
    CGAL_precondition(!empty());
    return tree().bbox();
  }

  template <typename T>
  bool do_intersect(const T& t) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::do_intersect(t, tree());
  }

  FT squared_distance(const Point_3& p) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::squared_distance(p, tree());
  }

  Point_3 closest_point(const Point_3& p) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::closest_point(p, tree());
  }

  bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o,
                          const FT offset_size, const FT intersection_precision) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
  }

  bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o, const FT offset_size) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::first_intersection(p, q, o, offset_size, 1e-2 * offset_size, tree());
  }

private:
  Geom_traits m_gt;
  AABB_tree_ptr m_tree_ptr;
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H
