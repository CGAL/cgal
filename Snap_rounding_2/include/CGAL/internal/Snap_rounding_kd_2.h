// Copyright (c) 2001, 2009, 2014 Tel-Aviv University (Israel), Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s) : Eli Packer <elip@post.tau.ac.il>,
//             Waqar Khan <wkhan@mpi-inf.mpg.de>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SNAP_ROUNDING_KD_2_H
#define CGAL_SNAP_ROUNDING_KD_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Dimension.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Snap_Rounding_2::internal {

/*! the interval number type; its operations require the rounding mode to be upward. */
using Interval = Interval_nt<false>;

/*! a point of the kd-tree, which represents a hot-pixel center.
 */
template <typename SavedObject>
struct Hot_pixel_kd_point {
  std::array<double, 2> coords;         // the midpoints of `x` and `y`, used by the kd-tree
  Interval x, y;                        // interval approximations of the center coordinates
  SavedObject object;                   // the hot pixel
};

/*! the search traits of the kd-tree of the hot pixels, a model of `SearchTraits`.
 */
template <typename SavedObject>
struct Hot_pixel_kd_search_traits {
  using Dimension = Dimension_tag<2>;
  using FT = double;
  using Point_d = Hot_pixel_kd_point<SavedObject>;
  using Cartesian_const_iterator_d = const double*;

  struct Construct_cartesian_const_iterator_d {
    Cartesian_const_iterator_d operator()(const Point_d& p) const { return p.coords.data(); }
    Cartesian_const_iterator_d operator()(const Point_d& p, int) const { return p.coords.data() + 2; }
  };

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const { return {}; }
};

/*! a fuzzy query item, a model of `FuzzyQueryItem`, that contains the Minkowski sum of a segment
 * and a pixel. The Minkowski sum is a hexagon; it is the intersection of the bounding box
 * of the segment enlarged by the half pixel size \f$h\f$, and the slab
 * \f$|a(x-s_x) + b(y-s_y)| \le h(|a|+|b|)\f$, where \f$(a,b)\f$ is the normal of the
 * segment. The tests use interval arithmetic and reject only what is certainly
 * outside, so the query reports a superset of the hot pixels the segment
 * intersects. A point may lie outside the kd-tree cell that contains the
 * midpoint of its intervals, but by at most `point_error`, so the cells are
 * enlarged by `point_error`. The rounding mode must be upward during the query.
 */
template <typename SearchTraits>
class Segment_pixel_corridor {
public:
  using Dimension = typename SearchTraits::Dimension;
  using FT = double;
  using Point_d = typename SearchTraits::Point_d;
  using Rectangle = Kd_tree_rectangle<FT, Dimension>;

private:
  Interval m_sx, m_sy;          // the source
  Interval m_a, m_b;            // the normal of the segment
  double m_w;                   // an upper bound on h*(|a|+|b|)
  Interval m_point_error;
  std::array<double, 2> m_min;  // the enlarged bounding box
  std::array<double, 2> m_max;

  /*! determines whether all the points in \f$x \times y\f$ are outside the slab.
   */
  bool certainly_outside_slab(const Interval& x, const Interval& y) const {
    Interval f = m_a * (x - m_sx) + m_b * (y - m_sy);
    return (f.inf() > m_w) || (f.sup() < -m_w);
  }

public:
  /*! constructs the corridor of a segment.
   * \param sx, sy the coordinates of the source.
   * \param tx, ty the coordinates of the target.
   * \param pixel_size the pixel size.
   * \param point_error the maximum width of the intervals of the points in the kd-tree.
   */
  Segment_pixel_corridor(const Interval& sx, const Interval& sy, const Interval& tx, const Interval& ty,
                         const Interval& pixel_size, double point_error) :
    m_sx(sx), m_sy(sy), m_a(sy - ty), m_b(tx - sx), m_point_error(point_error) {
    const Interval h = pixel_size / 2;
    m_w = (h * (CGAL::abs(m_a) + CGAL::abs(m_b))).sup();
    m_min = {(Interval((std::min)(sx.inf(), tx.inf())) - h).inf(),
             (Interval((std::min)(sy.inf(), ty.inf())) - h).inf()};
    m_max = {(Interval((std::max)(sx.sup(), tx.sup())) + h).sup(),
             (Interval((std::max)(sy.sup(), ty.sup())) + h).sup()};
  }

  /*! determines whether a point may be in the corridor.
   */
  bool contains(const Point_d& p) const {
    if ((p.x.sup() < m_min[0]) || (p.x.inf() > m_max[0]) || (p.y.sup() < m_min[1]) || (p.y.inf() > m_max[1]))
      return false;
    return ! certainly_outside_slab(p.x, p.y);
  }

  /*! determines whether a kd-tree cell may contain a point in the corridor.
   */
  bool inner_range_intersects(const Rectangle& rect) const {
    const double xmin = (std::max)((rect.min_coord(0) - m_point_error).inf(), m_min[0]);
    const double xmax = (std::min)((rect.max_coord(0) + m_point_error).sup(), m_max[0]);
    const double ymin = (std::max)((rect.min_coord(1) - m_point_error).inf(), m_min[1]);
    const double ymax = (std::min)((rect.max_coord(1) + m_point_error).sup(), m_max[1]);
    if ((xmin > xmax) || (ymin > ymax)) return false;
    return ! certainly_outside_slab(Interval(xmin, xmax), Interval(ymin, ymax));
  }

  /*! determines whether a kd-tree cell is contained in the corridor. A wrong positive answer only
   * reports extra candidates.
   */
  bool outer_range_contains(const Rectangle& rect) const {
    if ((rect.min_coord(0) < m_min[0]) || (rect.max_coord(0) > m_max[0]) ||
        (rect.min_coord(1) < m_min[1]) || (rect.max_coord(1) > m_max[1]))
      return false;
    Interval f = m_a * (Interval(rect.min_coord(0), rect.max_coord(0)) - m_sx) +
      m_b * (Interval(rect.min_coord(1), rect.max_coord(1)) - m_sy);
    return (f.inf() >= -m_w) && (f.sup() <= m_w);
  }
};

/*! a kd-tree of the hot pixels. It does not own the saved objects.
 * The tree stores double approximations of the hot-pixel centers; a query
 * returns a superset of the hot pixels a segment intersects. The exact test
 * is the responsibility of the caller.
 */
template <typename Traits_, typename SavedObject>
class Hot_pixel_kd_tree {
  static_assert(std::is_pointer<SavedObject>::value, "SavedObject is not a pointer.");

private:
  using Traits = Traits_;
  using NT = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;

  using Search_traits = Hot_pixel_kd_search_traits<SavedObject>;
  using Kd_point = typename Search_traits::Point_d;
  using Kd_tree = CGAL::Kd_tree<Search_traits>;
  using Corridor = Segment_pixel_corridor<Search_traits>;

  struct Object_inserter {
    std::vector<SavedObject>* result;
    void operator()(const Kd_point& p) const { result->push_back(p.object); }
  };

  Kd_tree m_tree;
  Interval m_pixel_size;
  double m_point_error = 0;     // the maximum width of the point intervals

public:
  using Point_saved_pair_list = std::vector<std::pair<Point_2, SavedObject>>;

  /*! constructs the kd-tree.
   * \param points the hot-pixel centers and the corresponding objects.
   * \param pixel_size the pixel size.
   */
  Hot_pixel_kd_tree(const Point_saved_pair_list& points, NT pixel_size) :
    m_pixel_size(CGAL::to_interval(pixel_size)) {
    m_tree.reserve(points.size());
    for (const auto& [p, obj] : points) {
      Interval x(CGAL::to_interval(p.x())), y(CGAL::to_interval(p.y()));
      m_point_error = (std::max)({m_point_error, x.sup() - x.inf(), y.sup() - y.inf()});
      m_tree.insert(Kd_point{{(x.inf() + x.sup()) / 2, (y.inf() + y.sup()) / 2}, x, y, obj});
    }
    m_tree.build();
  }

  Hot_pixel_kd_tree(const Hot_pixel_kd_tree&) = delete;
  Hot_pixel_kd_tree& operator=(const Hot_pixel_kd_tree&) = delete;

  /*! appends to `result` the objects of a superset of the hot pixels that `s` intersects.
   * \param result the output.
   * \param s the query segment.
   * \param gt the traits.
   */
  void get_intersecting_points(std::vector<SavedObject>& result, const Segment_2& s, const Traits& gt) const {
    auto ctr_vertex = gt.construct_vertex_2_object();
    const Point_2& src = ctr_vertex(s, 0);
    const Point_2& trg = ctr_vertex(s, 1);
    Protect_FPU_rounding<true> protection;
    Corridor corridor(CGAL::to_interval(src.x()), CGAL::to_interval(src.y()),
                      CGAL::to_interval(trg.x()), CGAL::to_interval(trg.y()), m_pixel_size, m_point_error);
    m_tree.search(boost::make_function_output_iterator(Object_inserter{&result}), corridor);
  }
};

} // namespace Snap_Rounding_2::internal
} // namespace CGAL

#endif
