// Copyright (c) 2026 INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Efi Fogel    <efifogel@gmail.com>

#ifndef CGAL_AABB_TREE_INTERNAL_AABB_CONSERVATIVE_BBOX_TEST_H
#define CGAL_AABB_TREE_INTERNAL_AABB_CONSERVATIVE_BBOX_TEST_H

#include <CGAL/license/AABB_tree.h>

/*! \file
 *
 * Conservative tests of bounding boxes against linear queries (lines, rays,
 * and segments), used to prune the traversal of an AABB tree.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

#include <CGAL/Bbox_3.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace internal {
namespace AABB_tree {

/*! \class
 *
 * Conservative_linear_bbox_test is a test of a bounding box against a linear
 * query (a line, a ray, or a segment), used to prune the traversal of an AABB
 * tree. The test is conservative: it may report an intersection with a box
 * that the query does not intersect, but it never misses a box that the query
 * intersects. The traversal applies the exact test of the traits to the
 * primitives, so the results of the queries do not change.
 *
 * The query is \f$p + t d\f$, where \f$t \in (-\infty, +\infty)\f$ (a
 * line), \f$t \in [0, +\infty)\f$ (a ray), or \f$t \in [0, 1]\f$ (a
 * segment). The coordinates of \f$p\f$ and \f$d\f$ are given as intervals
 * that contain the exact values (exact for a kernel with double coordinates;
 * the approximations for a lazy kernel), so that the test is conservative for
 * the exact query. The constructor computes the inverse directions once with
 * outward rounding. The test of a box computes, for every axis, the range of
 * the parameter in the slab of the box in plain floating-point arithmetic, and
 * inflates the result to cover the rounding errors:
 * - Every computed product has a relative error of at most about \f$3u\f$
 *   (\f$u = 2^{-53}\f$), in any rounding mode, unless an operation overflows
 *   or underflows. The relative inflation (\f$10^{-14}\f$) is much larger.
 * - An underflow causes an absolute error below \f$10^{-307}\f$, which the
 *   absolute inflation (\f$10^{-300}\f$) covers.
 * - A non-finite value (an overflow, or \f$0 \cdot \infty\f$) makes the
 *   test accept the box.
 */
class Conservative_linear_bbox_test {
public:
  /*! constructs the test of the query \f$p + t d\f$, where
   * \f$t \in [t_{\min}, t_{\max}]\f$.
   * \param p the intervals that contain the coordinates of \f$p\f$.
   * \param d the intervals that contain the coordinates of \f$d\f$.
   * \param bounded_0 indicates whether \f$t_{\min} = 0\f$; otherwise,
   *        \f$t_{\min} = -\infty\f$.
   * \param bounded_1 indicates whether \f$t_{\max} = 1\f$; otherwise,
   *        \f$t_{\max} = +\infty\f$.
   */
  Conservative_linear_bbox_test(const Interval_nt<>* p, const Interval_nt<>* d, bool bounded_0, bool bounded_1) :
    m_bounded_0(bounded_0), m_bounded_1(bounded_1) {
    typename Interval_nt<>::Protector protector;
    for (int i = 0; i < 3; ++i) {
      m_p_lo[i] = p[i].inf();
      m_p_hi[i] = p[i].sup();
      if (d[i].inf() == 0 && d[i].sup() == 0) m_kind[i] = FIXED;
      else if (d[i].inf() <= 0 && d[i].sup() >= 0) m_kind[i] = FREE;
      else {
        m_kind[i] = SLAB;
        Interval_nt<> inv = 1 / d[i];
        m_inv_lo[i] = inv.inf();
        m_inv_hi[i] = inv.sup();
        if (! std::isfinite(m_inv_lo[i]) || ! std::isfinite(m_inv_hi[i]))
          m_unsafe = true;
      }
      if (! std::isfinite(m_p_lo[i]) || ! std::isfinite(m_p_hi[i])) m_unsafe = true;
    }
  }

  /*! tests the bounding box against the query conservatively.
   * \param bb the bounding box.
   * \return `false` if the query does not intersect `bb`; `true` if it may
   *         intersect `bb` (always, if it does).
   */
  bool operator()(const Bbox_3& bb) const {
    if (m_unsafe) return true;
    double lo = m_bounded_0 ? 0. : -std::numeric_limits<double>::infinity();
    double hi = m_bounded_1 ? 1. : std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; ++i) {
      switch(m_kind[i]) {
       case FIXED:
        if (m_p_hi[i] < (bb.min)(i) || m_p_lo[i] > (bb.max)(i)) return false;
        break;

       case FREE: break;
       case SLAB:
        {
          // t = n * inv, where n = c - p for c in [bmin, bmax] and p in
          // [p_lo, p_hi]; the range of t is spanned by the products of the
          // bounds of n and inv.
          const double n_lo = (bb.min)(i) - m_p_hi[i];
          const double n_hi = (bb.max)(i) - m_p_lo[i];
          const double a = n_lo * m_inv_lo[i];
          const double b = n_lo * m_inv_hi[i];
          const double c = n_hi * m_inv_lo[i];
          const double d = n_hi * m_inv_hi[i];
          if (! std::isfinite(a) || ! std::isfinite(b) || ! std::isfinite(c) || ! std::isfinite(d)) return true;
          lo = (std::max)(lo, (std::min)((std::min)(a, b), (std::min)(c, d)));
          hi = (std::min)(hi, (std::max)((std::max)(a, b), (std::max)(c, d)));
        }
        break;
      }
    }
    constexpr double eps = 1e-14;
    constexpr double tiny = 1e-300;
    return (lo - eps * std::abs(lo) - tiny) <= (hi + eps * std::abs(hi) + tiny);
  }

private:
  /*! the kind of an axis: the coordinate of \f$d\f$ is 0 (FIXED), its
   * interval contains 0 (FREE), or it does not (SLAB).
   */
  enum Kind { FIXED, FREE, SLAB };

  double m_p_lo[3], m_p_hi[3];                                          //!< the intervals of \f$p\f$
  double m_inv_lo[3] = { 0., 0., 0. }, m_inv_hi[3] = { 0., 0., 0. };  //!< the intervals of \f$1/d\f$
  Kind m_kind[3];                                                       //!< the kinds of the axes
  bool m_bounded_0;                                                     //!< is \f$t_{\min} = 0\f$?
  bool m_bounded_1;                                                     //!< is \f$t_{\max} = 1\f$?
  bool m_unsafe = false;                                                //!< accept every box?
};

/*! \class
 *
 * Bbox_pruning_test is the test of a bounding box against a query used to prune
 * the traversal of an AABB tree. By default, it is the exact test of the
 * traits.
 * \tparam AABBTraits the traits of the AABB tree.
 * \tparam Query the type of the query.
 */
template <typename AABBTraits, typename Query, typename = void>
class Bbox_pruning_test {
public:
  /*! constructs.
   * \param traits the traits of the AABB tree.
   * \param query the query; it must outlive the test.
   */
  Bbox_pruning_test(const AABBTraits& traits, const Query& query) : m_traits(traits), m_query(query) {}

  /*! tests the bounding box against the query.
   * \param bb the bounding box.
   * \return `true` if the query intersects `bb`; `false` otherwise.
   */
  bool operator()(const typename AABBTraits::Bounding_box& bb) const {
    return m_traits.do_intersect_object()(m_query, bb);
  }

private:
  const AABBTraits& m_traits;
  const Query& m_query;
};

/*! converts a number to an interval that contains it.
 */
template <typename T>
Interval_nt<> to_interval_nt(const T& x) { return Interval_nt<>(CGAL::to_interval(x)); }

/*! \class
 *
 * The specialization of Bbox_pruning_test for a line, a ray, or a segment of
 * the geometric traits of AABB traits whose bounding box is `Bbox_3`: the
 * conservative test (see Conservative_linear_bbox_test).
 */
template <typename AABBTraits, typename Query>
class Bbox_pruning_test<AABBTraits, Query,
                        std::enable_if_t<std::is_same_v<typename AABBTraits::Bounding_box, Bbox_3> &&
                                         (std::is_same_v<Query, typename AABBTraits::Geom_traits::Line_3> ||
                                          std::is_same_v<Query, typename AABBTraits::Geom_traits::Ray_3> ||
                                          std::is_same_v<Query, typename AABBTraits::Geom_traits::Segment_3>)>> {
  using GT = typename AABBTraits::Geom_traits;

public:
  /*! constructs; prepares the conservative test of the query.
   * \param query the query.
   */
  Bbox_pruning_test(const AABBTraits&, const Query& query) : m_test(make(query)) {}

  /*! tests the bounding box against the query conservatively.
   * \param bb the bounding box.
   * \return `false` if the query does not intersect `bb`; `true` if it may
   *         intersect `bb` (always, if it does).
   */
  bool operator()(const Bbox_3& bb) const { return m_test(bb); }

private:
  /*! sets the intervals that contain the coordinates of a point or a vector.
   */
  template <typename Point>
  static void set(Interval_nt<>* x, const Point& p) {
    x[0] = to_interval_nt(p.x());
    x[1] = to_interval_nt(p.y());
    x[2] = to_interval_nt(p.z());
  }

  /*! makes the conservative test of the query: a line \f$p + t v\f$, where
   * \f$p\f$ and \f$v\f$ are the point and the vector that the line stores,
   * or a ray or a segment \f$p + t (q - p)\f$, where \f$p\f$ and \f$q\f$
   * are its defining points.
   */
  static Conservative_linear_bbox_test make(const Query& query) {
    Interval_nt<> p[3], q[3], d[3];
    bool bounded_0 = false, bounded_1 = false;
    if constexpr(std::is_same_v<Query, typename GT::Line_3>) {
      // The point and the vector that the line stores.
      set(p, query.point());
      set(d, query.to_vector());
      return Conservative_linear_bbox_test(p, d, false, false);
    }
    else {
      if constexpr(std::is_same_v<Query, typename GT::Ray_3>) {
        set(p, query.source());
        set(q, query.second_point());
        bounded_0 = true;
      }
      else {
        set(p, query.source());
        set(q, query.target());
        bounded_0 = bounded_1 = true;
      }
      typename Interval_nt<>::Protector protector;
      for (int i = 0; i < 3; ++i) d[i] = q[i] - p[i];
    }
    return Conservative_linear_bbox_test(p, d, bounded_0, bounded_1);
  }

  Conservative_linear_bbox_test m_test;  //!< the conservative test
};

/*! \class
 *
 * Lazy_bbox_pruning_test is the pruning test of a traversal of an AABB tree. It
 * is prepared from the query at the first test of a node (the query does not
 * change during a traversal).
 * \tparam AABBTraits the traits of the AABB tree.
 * \tparam Query the type of the query.
 */
template <typename AABBTraits, typename Query>
class Lazy_bbox_pruning_test {
public:
  /*! tests the bounding box of a node against the query; prepares the test at
   * the first call.
   * \param traits the traits of the AABB tree.
   * \param query the query; the same in all calls.
   * \param bb the bounding box.
   * \return `false` if the query does not intersect `bb`; `true` if it may
   *         intersect `bb` (always, if it does).
   */
  bool operator()(const AABBTraits& traits, const Query& query, const typename AABBTraits::Bounding_box& bb) const {
    if (! m_test) m_test.emplace(traits, query);
    return (*m_test)(bb);
  }

private:
  mutable std::optional<Bbox_pruning_test<AABBTraits, Query>> m_test;  //!< the test, once prepared
};

} // namespace AABB_tree
} // namespace internal
} // namespace CGAL

#endif
