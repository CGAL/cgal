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

#ifndef CGAL_CANVAS_BASE_CANVAS_H
#define CGAL_CANVAS_BASE_CANVAS_H

#include <CGAL/Canvas/Base_canvas_point.h>
#include <CGAL/Canvas/Canvas_seeds.h>
#include <CGAL/Canvas/Metric.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <Eigen/Dense>
#error
#endif

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <vector>
#include <utility>

namespace CGAL {
namespace Canvas {

template<typename Cp>
struct Canvas_point_comparer
{
  bool operator()(Cp const * const cp1, Cp const * const cp2)
  {
    return cp1->distance_to_closest_seed() > cp2->distance_to_closest_seed();
  }
};

template<typename GeomTraits, typename CanvasPoint, typename MetricField>
class Base_canvas
{
private:
  using Self = Base_canvas<GeomTraits, CanvasPoint, MetricField>;

public:
  using Geom_traits = GeomTraits;

  using Metric_field = MetricField;
  using Metric = Metric_base<Geom_traits>;

  using Canvas_point = CanvasPoint;
  using Canvas_point_handle_vector = std::vector<Canvas_point*>;
  using Seeds = Canvas_seeds<Self>;

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector3d = Eigen::Matrix<FT, 3, 1>;

protected:
  std::vector<Canvas_point> m_canvas_points;
  Canvas_point_handle_vector m_trial_points;
  Seeds m_seeds;
  const Metric_field* m_mf;
  Geom_traits m_gt;

  // debug & info
  std::size_t m_known_count, m_trial_count, m_far_count;

public:
  Base_canvas(const Metric_field* mf_,
              const Geom_traits& gt = Geom_traits())
    :
      m_canvas_points(),
      m_trial_points(),
      m_seeds(*this),
      m_mf(mf_),
      m_known_count(0),
      m_trial_count(0),
      m_far_count(0),
      m_gt(gt)
  {
  }

public:
  const Metric_field* metric_field() const { return m_mf; }
  const std::vector<Canvas_point>& canvas_points() const { return m_canvas_points; }

  virtual void initialize() = 0;

  Canvas_point& get_canvas_point(const std::size_t i)
  {
    CGAL_precondition(i < m_canvas_points.size());
    return m_canvas_points[i];
  }

  const Canvas_point& get_canvas_point(const std::size_t i) const
  {
    CGAL_precondition(i < m_canvas_points.size());
    return m_canvas_points[i];
  }

  virtual bool is_point_outside_canvas(const FT, const FT, const FT) const
  {
    return false;
  }

  virtual bool is_point_outside_canvas(const Point_3& p) const
  {
    return is_point_outside_canvas(p.x(), p.y(), p.z());
  }

  void initialize_canvas_point(Canvas_point& cp,
                               const FT distance_from_seed,
                               const std::size_t seed_id)
  {
    if(cp.closest_seed_id() != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed is overwriting the closest seed id";
      std::cout << " of a canvas point!" << std::endl;
      std::cout << "seed: " << seed_id << " wants to initialize point: " << cp.index() << "(" << cp.point() << ")";
      std::cout << " but seed " << cp.closest_seed_id() << " has done it already" << std::endl;
    }

    // We can't accept two seeds for one canvas point
    if(cp.state() == TRIAL)
      CGAL_assertion(false && "the canvas is not dense enough for the input seeds...");

    cp.initialize_from_point(distance_from_seed, seed_id);

    m_trial_points.push_back(&cp);
    std::push_heap(m_trial_points.begin(), m_trial_points.end(),
                   Canvas_point_comparer<Canvas_point>());
  }

  virtual void locate_and_initialize(const Point_3& s,
                                     const std::size_t seed_id)
  {
    // something pretty todo...
    // for now: brutally find the canvas point closest to the seed

    int cp_id = -1; // index of the closest seed

    FT min_d = std::numeric_limits<double>::infinity();
    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
    {
      Vector3d v;
      v(0) = s.x() - m_canvas_points[i].point().x();
      v(1) = s.y() - m_canvas_points[i].point().y();
      v(2) = s.z() - m_canvas_points[i].point().z();
      const Eigen::Matrix3d& m = m_canvas_points[i].metric().get_mat();

      FT d = v.transpose() * m * v; // note that this is the squared_distance
      if(d < min_d)
      {
        min_d = d;
        cp_id = i;
      }
    }

    CGAL_assertion(cp_id >= 0);
    CGAL_assertion(min_d == 0); // @tmp only using seeds that are vertices right now
    CGAL_assertion(min_d == CGAL::squared_distance(s, m_canvas_points[cp_id].point())); // @tmp using Euclidean distance currently

    Canvas_point& cp = m_canvas_points[cp_id];

#if (VERBOSITY > 15)
    std::cout << "looking for seed: " << s << std::endl;
    std::cout << "found cp: " << cp.index() << " [" << cp.point() << "] ";
    std::cout << "at distance: " << std::sqrt(min_d) << std::endl;
#endif

    initialize_canvas_point(cp, std::sqrt(min_d), seed_id);
  }

  void locate_seeds_on_canvas()
  {
    // find the canvas vertex closest to the seed
    for(std::size_t i=0; i<m_seeds.size(); ++i)
    {
      const Point_3& p = m_seeds[i];
      locate_and_initialize(p, i);
    }
  }

  void print_states() const
  {
    std::cout << "known: " << m_known_count;
    std::cout << " trial: " << m_trial_count;
    std::cout << " far: " << m_far_count << std::endl;
  }

  void reset_counters()
  {
    m_known_count = 0;
    m_trial_count = 0;
    m_far_count = m_canvas_points.size();
  }

  void refresh_canvas_point_states()
  {
    CGAL_assertion(m_trial_points.empty());
    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
      m_canvas_points[i].state() = FAR;

    reset_counters();
  }

  void set_points_states_to_known()
  {
    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
      m_canvas_points[i].state() = KNOWN;

    m_known_count = m_canvas_points.size();
    m_trial_count = 0;
    m_far_count = 0;
  }

  void reset()
  {
    // Remove everything related to "paint" but not the geometric information
    // (therefore the point position, metric, neighbors, etc. aren't reset)

    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
    {
      Canvas_point& cp = m_canvas_points[i];
      cp.state() = FAR;
      cp.distance_to_closest_seed() = std::numeric_limits<double>::infinity();
      cp.closest_seed_id() = -1;
      cp.ancestor() = -1;
      cp.children().clear();
    }

    reset_counters();
  }

  virtual void paint()
  {
#if (VERBOSITY > 5)
    std::cout << "Paiting..." << std::endl;
    std::clock_t start = std::clock();

    if(m_trial_points.empty())
      std::cerr << "Trying to paint without anything in the PQ..." << std::endl;
    else
      std::cout << m_trial_points.size() << " initial points in the queue" << std::endl;
#endif

    while(!m_trial_points.empty())
    {
#if (VERBOSITY > 10)
      if(known_count % ((std::max)(static_cast<std::size_t>(1), m_canvas_points.size()/100)) == 0)
        print_states();
#endif

#if (VERBOSITY > 15)
      std::cout << "Trial queue size : " << trial_points.size() << std::endl;
#endif

#if (VERBOSITY > 55)
      std::cout << "trial heap: " << std::endl;
      for(Canvas_point* tp : m_trial_points)
        std::cout << tp->index() << " " << tp->distance_to_closest_seed() << std::endl;
      std::cout << std::endl;
#endif

      Canvas_point* cp = m_trial_points.front();
      CGAL_assertion(cp);
      CGAL_assertion(cp->state() == TRIAL);
      std::pop_heap(m_trial_points.begin(), m_trial_points.end(), Canvas_point_comparer<Canvas_point>());
      m_trial_points.pop_back();

#if (VERBOSITY > 15)
      std::cout << "picked n° " << cp->index() << " (" << cp->point() << ") ";
      std::cout << "at distance : " << cp->distance_to_closest_seed()
                << " from " << cp->closest_seed_id() << std::endl;
#endif

      cp->change_state(KNOWN, m_known_count, m_trial_count, m_far_count);
      PQ_state pqs = cp->update_neighbors_distances(m_trial_points);

      if(pqs == REBUILD_TRIAL)
        std::make_heap(m_trial_points.begin(), m_trial_points.end(), Canvas_point_comparer<Canvas_point>());
    }

    CGAL_expensive_assertion_code(debug());

#if (VERBOSITY > 15)
    std::cout << "final states after painting: " << std::endl;
    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
    {
      const Canvas_point& cp = m_canvas_points[i];
      std::cout << cp.index() << " (" << cp.point() << ")";
      std::cout << " at distance: " << cp.distance_to_closest_seed();
      std::cout << " from " << cp.closest_seed_id() << std::endl;
//      cp.print_ancestor_tree();
    }
#endif

#if (VERBOSITY > 5)
    std::cout << "End of paint. time: ";
    std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
#endif
  }

  virtual void debug()
  {
    CGAL_assertion(m_trial_points.empty() );

    for(std::size_t i=0, cps=m_canvas_points.size(); i<cps; ++i)
    {
      const Canvas_point& cp = m_canvas_points[i];

      if(cp.ancestor() != static_cast<std::size_t>(-1) &&
         m_canvas_points[cp.ancestor()].children().find(cp.index()) ==
         m_canvas_points[cp.ancestor()].children().end())
      {
        std::cout << "failure in ancestor/children relationship at " << cp.index() << std::endl;
        CGAL_assertion(false);
      }

      CGAL_postcondition(cp.distance_to_closest_seed() != std::numeric_limits<double>::infinity());
      CGAL_postcondition(cp.closest_seed_id() < m_seeds.size());
    }
  }
};

} // namespace Canvas
} // namespace CGAL

#endif // CGAL_CANVAS_BASE_CANVAS_H
