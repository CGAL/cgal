// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)
//                 Clement Jamin (clement.jamin.pro@gmail.com)


#ifndef CGAL_EUCLIDEAN_DISTANCE_H
#define CGAL_EUCLIDEAN_DISTANCE_H

#include <CGAL/license/Spatial_searching.h>


#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/number_utils.h>
#include <CGAL/internal/Get_dimension_tag.h>
#include <vector>
#include <iterator>

namespace CGAL {

  template <class SearchTraits>
  class Euclidean_distance;


  namespace internal{
    template <class SearchTraits>
    struct Spatial_searching_default_distance{
      typedef ::CGAL::Euclidean_distance<SearchTraits> type;
    };
  } //namespace internal

  template <class SearchTraits>
  class Euclidean_distance {

    SearchTraits traits;

    public:

    typedef typename SearchTraits::FT    FT;
    typedef typename SearchTraits::Point_d Point_d;
    typedef Point_d Query_item;

    typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension D;


    // default constructor
    Euclidean_distance(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}

    inline FT transformed_distance(const Query_item& q, const Point_d& p) const
    {
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d p_begin = construct_it(p), p_end = construct_it(p, 0);
      return transformed_distance_from_coordinates(q, p_begin, p_end);
    }

    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator it_coord_end) const
    {
      return transformed_distance_from_coordinates(q, it_coord_begin, it_coord_end, D());
    }

    // Static dim = 2 loop unrolled
    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
                                                    Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
                                                    Dimension_tag<2>) const
    {
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q);
      FT distance = square(*qit - *it_coord_begin);
      qit++; it_coord_begin++;
      distance += square(*qit - *it_coord_begin);
      return distance;
    }

    // Static dim = 3 loop unrolled
    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
                                                    Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
                                                    Dimension_tag<3>) const
    {
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q);
      FT distance = square(*qit - *it_coord_begin);
      qit++; it_coord_begin++;
      distance += square(*qit - *it_coord_begin);
      qit++; it_coord_begin++;
      distance += square(*qit - *it_coord_begin);
      return distance;
    }

    // Other cases: static dim > 3 or dynamic dim
    template <typename Coord_iterator, typename Dim>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
                                                    Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
                                                    Dim) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q), qe = construct_it(q, 1);
      for (; qit != qe; ++qit, ++it_coord_begin)
      {
        FT diff = (*qit) - (*it_coord_begin);
        distance += diff*diff;
      }
      return distance;
    }

    // During the computation, if the partially-computed distance `pcd` gets greater or equal
    // to `stop_if_geq_to_this`, the computation is stopped and `pcd` is returned
    template <typename Coord_iterator>
    inline FT interruptible_transformed_distance(const Query_item& q,
                                                 Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
                                                 FT stop_if_geq_to_this) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q), qe = construct_it(q, 1);
      if (std::distance(qit, qe) >= 6)
      {
        // Every 4 coordinates, the current partially-computed distance
        // is compared to stop_if_geq_to_this
        // Note: the concept SearchTraits specifies that Cartesian_const_iterator_d
        //       must be a random-access iterator
        typename SearchTraits::Cartesian_const_iterator_d qe_minus_5 = qe;
        std::advance(qe_minus_5, -5);
        for (;;)
        {
          FT diff = (*qit) - (*it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = (*qit) - (*it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = (*qit) - (*it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = (*qit) - (*it_coord_begin);
          distance += diff*diff;
          ++qit, ++it_coord_begin;

          if (distance >= stop_if_geq_to_this)
            return distance;

          if (std::distance(qe_minus_5, qit) >= 0)
            break;
        }
      }
      for (; qit != qe; ++qit, ++it_coord_begin)
      {
        FT diff = (*qit) - (*it_coord_begin);
        distance += diff*diff;
      }
      return distance;
    }

    inline FT min_distance_to_rectangle(const Query_item& q,
      const Kd_tree_rectangle<FT, D>& r) const {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
        qe = construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        if ((*qit) < r.min_coord(i))
          distance +=
          (r.min_coord(i) - (*qit))*(r.min_coord(i) - (*qit));
        else if ((*qit) > r.max_coord(i))
          distance +=
          ((*qit) - r.max_coord(i))*((*qit) - r.max_coord(i));

      }
      return distance;
    }

    inline FT min_distance_to_rectangle(const Query_item& q,
      const Kd_tree_rectangle<FT, D>& r, std::vector<FT>& dists) const {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
        qe = construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        if ((*qit) < r.min_coord(i)) {
          dists[i] = (r.min_coord(i) - (*qit));
          distance += dists[i] * dists[i];
        }
        else if ((*qit) > r.max_coord(i)) {
          dists[i] = ((*qit) - r.max_coord(i));
          distance += dists[i] * dists[i];
        }

      }
      return distance;
    }

    inline FT max_distance_to_rectangle(const Query_item& q,
      const Kd_tree_rectangle<FT, D>& r) const {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
        qe = construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        if ((*qit) <= (r.min_coord(i) + r.max_coord(i)) / FT(2.0))
          distance += (r.max_coord(i) - (*qit))*(r.max_coord(i) - (*qit));
        else
          distance += ((*qit) - r.min_coord(i))*((*qit) - r.min_coord(i));
      };
      return distance;
    }

    inline FT max_distance_to_rectangle(const Query_item& q,
      const Kd_tree_rectangle<FT, D>& r, std::vector<FT>& dists) const {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Cartesian_const_iterator_d qit = construct_it(q),
        qe = construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        if ((*qit) <= (r.min_coord(i) + r.max_coord(i)) / FT(2.0)) {
          dists[i] = (r.max_coord(i) - (*qit));
          distance += dists[i] * dists[i];
        }
        else {
          dists[i] = ((*qit) - r.min_coord(i));
          distance += dists[i] * dists[i];
        }
      };
      return distance;
    }

    inline FT new_distance(FT dist, FT old_off, FT new_off,
      int /* cutting_dimension */)  const {

      FT new_dist = dist + (new_off*new_off - old_off*old_off);
      return new_dist;
    }

    inline FT transformed_distance(FT d) const {
      return d*d;
    }

    inline FT inverse_of_transformed_distance(FT d) const {
      return CGAL::sqrt(d);
    }

  }; // class Euclidean_distance

} // namespace CGAL
#endif // EUCLIDEAN_DISTANCE_H
