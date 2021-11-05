// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_CANVAS_BASE_CANVAS_POINT_H
#define CGAL_CANVAS_BASE_CANVAS_POINT_H

#include <CGAL/Canvas/Metric.h>

#include <CGAL/assertions.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <Eigen/Dense>
#else
#error
#endif

#include <boost/unordered_set.hpp>

#include <cstddef>
#include <iostream>

namespace CGAL {
namespace Canvas {

enum FMM_state
{
  CHANGED = 0, // used only in the Konukoglu algorithm. 'CHANGED' implies 'KNOWN'
  KNOWN,
  TRIAL,
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL,
  REBUILD_CHANGED, // used only in the Konukoglu algorithm
  REBUILD_BOTH // used only in the Konukoglu algorithm
};

template<typename GeomTraits, typename Canvas>
class Base_canvas_point
{
private:
  using Self = Base_canvas_point<GeomTraits, Canvas>;

protected:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Metric = Metric_base<Geom_traits>;
  using Vector3d = Eigen::Matrix<FT, 3, 1>;

  using Point_set = boost::unordered_set<std::size_t>;

protected:
  const Point_3& m_point;
  Canvas& m_canvas;

  std::size_t m_index;
  Metric m_metric;

  // stuff that depends on the seeds :
  FMM_state m_state;
  FT m_distance_to_closest_seed;
  std::size_t m_depth;
  int m_is_seed_holder; // ugly and awkward (should use the fact that it has no ancestor instead)
  std::size_t m_closest_seed_id;
  std::size_t m_ancestor;

  // 'children' needs to be 'mutable' because compute_closest_seed takes a const ref
  // to an ancestor (because some ancestors are temporaries whose lifestime I extend
  // through const refs) and we need to modify children in compute_closest_seed
  mutable Point_set m_children;

public:
  //  Base_canvas_point() { }
  Base_canvas_point(const std::size_t index_,
                    const Point_3& point_,
                    Canvas& canvas)
      :
        m_canvas(canvas),
        m_point(point_),
        m_index(index_),
        m_metric(canvas.metric_field()->compute_metric(m_point)),
        m_state(FAR),
        m_distance_to_closest_seed(std::numeric_limits<double>::infinity()),
        m_depth(0),
        m_is_seed_holder(-1),
        m_closest_seed_id(-1),
        m_ancestor(-1),
        m_children()
    { }

public:
  const Point_3& point() const { return m_point; }
  std::size_t& index() { return m_index; }
  std::size_t index() const { return m_index; }
  FT& distance_to_closest_seed() { return m_distance_to_closest_seed; }
  const FT distance_to_closest_seed() const { return m_distance_to_closest_seed; }
  int& is_seed_holder() { return m_is_seed_holder; }
  int is_seed_holder() const { return m_is_seed_holder; }
  std::size_t& depth() { return m_depth; }
  std::size_t depth() const { return m_depth; }
  std::size_t& closest_seed_id() { return m_closest_seed_id; }
  std::size_t closest_seed_id() const { return m_closest_seed_id; }
  FMM_state& state() { return m_state; }
  FMM_state state() const { return m_state; }
  Metric& metric() { return m_metric; }
  const Metric& metric() const { return m_metric; }
  std::size_t& ancestor() { return m_ancestor; }
  std::size_t ancestor() const { return m_ancestor; }
  Canvas& canvas() { return m_canvas; }
  const Canvas& canvas() const { return m_canvas; }
  Point_set& children() { return m_children; }
  const Point_set& children() const { return m_children; }

  void change_state(FMM_state new_state, std::size_t& known_count,
                    std::size_t& trial_count, std::size_t& far_count)
  {
    if(new_state == m_state)
      std::cerr << "WARNING: useless state change..." << std::endl;

    if(state() == FAR)
      --far_count;
    else if(state() == TRIAL)
      --trial_count;
    else if(state() == KNOWN)
      --known_count;

    if(new_state == KNOWN)
      ++known_count;
    else if(new_state == TRIAL)
      ++trial_count;
    else if(new_state == FAR)
      ++far_count;

    state() = new_state;
  }

  void remove_from_children(const std::size_t c)
  {
    // 99% of the time, the child is found, but if during the initialization of
    // a new seed, you reset (at least) two points that have an ancestry relationship
    // then the child could have been reset already...
    typename Point_set::iterator it = m_children.find(c);
    if(it != m_children.end())
      m_children.quick_erase(it);
#if (VERBOSITY > 15)
    else
    {
      std::cerr << "WARNING: call to remove_from_children didn't find the child (";
      std::cerr << c << " from " << index << ")" << std::endl;
    }
#endif
  }

  void reset_descendants()
  {
#if (VERBOSITY > 55)
     std::cout << "reset descendant at: " << index << std::endl;
#endif
    while(!m_children.empty())
    {
      Self& cp = canvas().get_canvas_point(*(m_children.begin()));
      cp.reset_descendants();
    }
    reset_color();
  }

  void initialize_from_point(const FT d,
                             const std::size_t seed_id)
  {
#if (VERBOSITY > 15)
    std::cout << "initialize " << m_index << " (" << m_point << ")";
    std::cout << " at distance " << d << " from " << seed_id << std::endl;
#endif
    reset_color();

    m_is_seed_holder = seed_id;
    m_closest_seed_id = seed_id;
    m_distance_to_closest_seed = d;
    m_state = TRIAL;
  }

  void print_ancestor_tree() const
  {
    std::cout << m_index << " has ancestors: ";
    std::size_t anc = m_ancestor;
    while(anc != static_cast<std::size_t>(-1))
    {
      std::cout << anc << " ";
      anc = canvas().get_canvas_point(anc).ancestor();
    }
    std::cout << std::endl;
  }

  FT distortion_to_seed() const
  {
    FT gamma = 1.;
//    std::cout << "init gamma: " << gamma << " " << index << std::endl;
    std::size_t curr = m_index;
    std::size_t anc = m_ancestor;

    while(anc != static_cast<std::size_t>(-1))
    {
      const Metric& m1 = canvas().get_canvas_point(anc).metric();
      const Metric& m2 = canvas().get_canvas_point(curr).metric();
      FT loc_gamma = m1.compute_distortion(m2);
//      std::cout << "loc: " << loc_gamma << " " << anc->index << std::endl;
#if 1
      gamma *= loc_gamma;
#else
      gamma = (std::max)(loc_gamma, gamma);
#endif
//      std::cout << "gamma: " << gamma << std::endl;
      anc = canvas().canvas_points[anc].ancestor();
    }

    gamma = (std::min)(25., gamma);

//    std::cout << "final gamma:" << gamma << std::endl;
    return gamma;
  }

  std::size_t count_ancestors() const
  {
    std::size_t i = 1;
    std::size_t anc = m_ancestor;
    while(anc != static_cast<std::size_t>(-1))
    {
      ++i;
      anc = canvas().get_canvas_point(anc).ancestor();
    }
    return i;
  }

  std::size_t ancestor_path_length() const
  {
    std::size_t i = 1;
    std::size_t n_anc = m_ancestor;
    while(n_anc != static_cast<std::size_t>(-1))
    {
      ++i;
      const Self& anc = canvas().get_canvas_point(n_anc);
      n_anc = anc.ancestor();

      CGAL_assertion(i < canvas().canvas_points().size());
    }
    return i;
  }

  void reset_color()
  {
    if(m_ancestor != static_cast<std::size_t>(-1))
      canvas().get_canvas_point(m_ancestor).remove_from_children(m_index);

    for(std::size_t cid : m_children)
      canvas().get_canvas_point(cid).ancestor() = -1;

    m_state = FAR;
    m_distance_to_closest_seed = std::numeric_limits<double>::infinity();
    m_depth = 0;
    m_is_seed_holder = -1;
    m_closest_seed_id = -1;
    m_ancestor = -1;
    m_children.clear();
  }
};

}  // namespace Canvas
}  // namespace CGAL

#endif // CGAL_CANVAS_BASE_CANVAS_POINT_H
