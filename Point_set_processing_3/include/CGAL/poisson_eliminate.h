// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau

#ifndef CGAL_POISSON_ELIMINATE_H
#define CGAL_POISSON_ELIMINATE_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>


namespace CGAL {

namespace internal {

double get_maximum_radius(std::size_t dimension, std::size_t sample_size, double domain_size) {
  double sampleArea = domain_size / double(sample_size);
  double r_max;
  switch (dimension) {
  case 2: r_max = CGAL::sqrt(sampleArea / (2.0 * CGAL::sqrt(3.0))); break;
  case 3: r_max = std::pow(sampleArea / (4.0 * CGAL::sqrt(2.0)), 1.0 / 3.0); break;
  default:
    double c;
    std::size_t d_start;
    if ((dimension & 1)) {
      c = 2.0;
      d_start = 3;
    }
    else {
      c = CGAL_PI;
      d_start = 4;
    }
    for (std::size_t d = d_start; d <= dimension; d += 2)
      c *= 2.0 * CGAL_PI / double(d);

    r_max = std::pow(sampleArea / c, 1.0 / double(dimension));
    break;
  }
  return r_max;
}

double get_minimum_radius(std::size_t input_size, std::size_t output_size, double beta, double gamma, double r_max) {
  double ratio = output_size / double(input_size);
  return r_max * (1 - std::pow(ratio, gamma)) * beta;
}

template<class Point, class FT = typename CGAL::Kernel_traits<Point>::FT, unsigned int dim = CGAL::Ambient_dimension<Point>::value>
Point construct(const std::array<FT, dim>& coords) {
  if constexpr (dim == 2) { return Point(coords[0], coords[1]); }
  else if constexpr (dim == 3) { return Point(coords[0], coords[1], coords[2]); }
  else { return Point(coords); }
}

template<class Point, class FT = typename CGAL::Kernel_traits<Point>::FT>
void copy_and_replace(const Point& in, Point& out, std::size_t replace, FT value) {
  std::array<FT, CGAL::Ambient_dimension<Point>::value> tmp;
  for (std::size_t i = 0; i < CGAL::Ambient_dimension<Point>::value; i++)
    if (i != replace)
      tmp[i] = in[int(i)];
    else
      tmp[replace] = value;

  out = construct<Point, FT, CGAL::Ambient_dimension<Point>::value>(tmp);
}

template<class PointRange, class PointMap>
class Indexed_extended_point_map {
public:
  using value_type = typename boost::property_traits<PointMap>::value_type;
  using reference = std::conditional_t< std::is_reference_v<typename boost::property_traits<PointMap>::reference>,
                                        const value_type&, value_type>;
  using key_type = std::size_t;
  using category = boost::readable_property_map_tag;

private:
  const PointRange& range;
  const std::vector<value_type>& tiling_points;
  const PointMap& point_map;

public:
  Indexed_extended_point_map(const PointRange& input, const std::vector<value_type>& tiling, const PointMap& point_map)
    : range(input), tiling_points(tiling), point_map(point_map) {}

  reference operator[](key_type k) const {
    CGAL_assertion(k < (range.size() + tiling_points.size()));
    if (range.size() > k)
      return get(point_map, *(range.begin() + k));
    else
      return tiling_points[k - range.size()];
  }

  friend reference get(const Indexed_extended_point_map& ppmap, key_type k) {
    CGAL_assertion(k < (ppmap.range.size() + ppmap.tiling_points.size()));
    if ((k < ppmap.range.size()))
      return get(ppmap.point_map, *(ppmap.range.begin() + k));
    else
      return ppmap.tiling_points[k - ppmap.range.size()];
  }
};

template<class Point>
class Weight_functor {
public:
  Weight_functor(double r_min = 0, double alpha = 8) : r_min(CGAL::to_double(r_min)), alpha(CGAL::to_double(alpha)) {}

  double operator()(const Point &, const Point &, double d2, double r_max) {
    double d = CGAL::sqrt(d2);
    if (d < r_min) d = r_min;
    return std::pow(double(1) - d / r_max, alpha);
  }

private:
  double r_min;
  double alpha;
};

void move_down(std::vector<std::size_t>& heap, std::vector<std::size_t>& heap_pos, std::size_t heap_size, std::vector<double>& weights, std::size_t idx) {
  CGAL_assertion(idx <= heap.size());
  std::size_t child = idx * 2 + 1;

  while (child + 1 < heap_size) {
    if (weights[heap[child]] < weights[heap[child + 1]]) child++;
    if (weights[heap[idx]] >= weights[heap[child]])
      break;

    std::swap(heap[idx], heap[child]);
    heap_pos[heap[idx]] = idx;
    heap_pos[heap[child]] = child;

    idx = child;
    child = idx * 2 + 1;
  }

  if (child < heap_size) {
    if (weights[heap[idx]] < weights[heap[child]]) {
      std::swap(heap[idx], heap[child]);

      heap_pos[heap[idx]] = idx;
      heap_pos[heap[child]] = child;
    }
  }
}

void pop_heap(std::vector<std::size_t>& heap, std::vector<std::size_t>& heap_pos, std::size_t &heap_size, std::vector<double>& weights) {
  std::swap(heap.front(), heap[heap_size - 1]);
  heap_pos[heap.front()] = 0;
  heap_pos[heap[heap_size - 1]] = heap_size - 1;

  heap_size--;

  move_down(heap, heap_pos, heap_size, weights, 0);
}

template<class Point, int d = CGAL::Ambient_dimension<Point>::value>
struct Compute_squared_distance {
  using GeomTraits = typename CGAL::Kernel_traits<Point>::Kernel;
  using FT = typename GeomTraits::FT;

  FT operator()(const Point& a, const Point& b) {
    return GeomTraits().squared_distance_d_object()(a, b);
  }
};

template<class Point>
struct Compute_squared_distance<Point, 2> {
  using GeomTraits = typename CGAL::Kernel_traits<Point>::Kernel;
  using FT = typename GeomTraits::FT;

  FT operator()(const Point& a, const Point& b) {
    return GeomTraits().compute_squared_distance_2_object()(a, b);
  }
};

template<class Point>
struct Compute_squared_distance<Point, 3> {
  using GeomTraits = typename CGAL::Kernel_traits<Point>::Kernel;
  using FT = typename GeomTraits::FT;

  FT operator()(const Point& a, const Point& b) {
    return GeomTraits().compute_squared_distance_3_object()(a, b);
  }
};

template<class Point, int d = CGAL::Ambient_dimension<Point>::value>
struct Search_traits : public CGAL::Search_traits_d<typename CGAL::Kernel_traits<Point>::Kernel>{
};

template<class Point>
struct Search_traits<Point, 2> : public CGAL::Search_traits_2<typename CGAL::Kernel_traits<Point>::Kernel> {
};

template<class Point>
struct Search_traits<Point, 3> : public CGAL::Search_traits_3<typename CGAL::Kernel_traits<Point>::Kernel> {
};

}

/**
   \ingroup PkgPointSetProcessing3Algorithms
   performs Poisson disk elimination with a desired output size. A greedy method that calculates a weight based on the
   neighborhood of each point and eliminates points until the output size is reached.

   For more details, please refer to \cgalCite{cgal:y-sefpdss}.

   \tparam PointRange is a model of `RandomAccessRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \tparam OutputIterator Type of the output iterator. Must accept input of the same type as the iterator of `PointRange`.

   \param points input point range
   \param number_of_points target output size, must be smaller that the number of points in the input range
   \param output where output points are put
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
       \cgalParamNBegin{point_map}
         \cgalParamDescription{a property map associating points to the elements of the point set `points`}
         \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                        of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
         \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
       \cgalParamNEnd

       \cgalParamNBegin{dimension}
         \cgalParamDescription{The sampling domain of `points`, e.g., 2 if the points have been sampled from a 2d surface.}
         \cgalParamType{unsigned integer}
         \cgalParamDefault{2}
       \cgalParamNEnd

       \cgalParamNBegin{progressive}
         \cgalParamDescription{reorders the points in `output` in a progressive way, i.e., the first n points in `output` with n < `number_of_points` have a Poisson disk distribution with a larger radius. }
         \cgalParamType{Boolean}
         \cgalParamDefault{`false`}
       \cgalParamNEnd

       \cgalParamNBegin{maximum_radius}
         \cgalParamDescription{radius of the Poisson disk in which the neighboring points are taken into account for elimination.}
         \cgalParamType{double}
         \cgalParamDefault{the default is calculated from the `dimension`, the volume of the bounding box and the `number_of_points`. For more details, see parameter section in \ref Point_set_processing_3PoissonElimination.}
       \cgalParamNEnd

       \cgalParamNBegin{weight_function}
         \cgalParamDescription{a weight function that calculates the weight of a neighbor point based on its squared distance and the `maximum_radius`.}
         \cgalParamType{an instance of `std::function<double(const Point&, const Point&, double,double)>`.}
         \cgalParamDefault{See parameter section in \ref Point_set_processing_3PoissonElimination.}
       \cgalParamNEnd

       \cgalParamNBegin{geom_traits}
         \cgalParamDescription{an instance of a geometric traits class}
         \cgalParamType{a model of `Kernel`}
         \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
       \cgalParamNEnd
     \cgalNamedParamsEnd
*/
template<class PointRange, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
void poisson_eliminate(const PointRange &points, std::size_t number_of_points, OutputIterator output, const NamedParameters& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  if (points.size() == 0)
    return;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters>;
  using PointMap = typename NP_helper::Const_point_map;
  using Point = typename boost::property_traits<PointMap>::value_type;
  using GeomTraits = typename NP_helper::Geom_traits;
  using FT = typename GeomTraits::FT;
  using IPM = internal::Indexed_extended_point_map<PointRange, PointMap>;
  PointMap point_map = NP_helper::get_const_point_map(points, np);

  const unsigned int ambient_dimension = CGAL::Ambient_dimension<Point>::value;

  using Search_base = internal::Search_traits<Point>;

  using Search_traits = CGAL::Search_traits_adapter<std::size_t, IPM, Search_base>;
  using Splitter = CGAL::Sliding_midpoint<Search_traits>;
  using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;
  using Tree = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true, CGAL::Tag_true>;

  internal::Compute_squared_distance<Point> squared_distance;

  // hard coded parameters alpha, beta and gamma with values proposed in publication
  const double alpha = 8;
  const double beta = 0.65; // not used currently
  const double gamma = 1.5; // not used currently
  // named parameters for weight limiting, currently not used as the performance gain seems small and there is a high risk to get a defective output
  const bool weight_limiting = parameters::choose_parameter(parameters::get_parameter(np, internal_np::weight_limiting), false);
  // named parameter for progressive
  const bool progressive = parameters::choose_parameter(parameters::get_parameter(np, internal_np::progressive), false);
  // named parameter for tiling
  const bool tiling = parameters::choose_parameter(parameters::get_parameter(np, internal_np::tiling), false);
  const unsigned int dimension = parameters::choose_parameter(parameters::get_parameter(np, internal_np::dimension), 2);

  CGAL_assertion(ambient_dimension >= dimension);

  // Needs a multi-dimension solution
  std::array<FT, ambient_dimension> lower, upper;
  const Point &first = get(point_map, *(points.begin()));
  for (unsigned int i = 0; i < ambient_dimension; i++)
    lower[i] = upper[i] = first[i];

  for (const auto &pt : points) {
    const Point& p = get(point_map, pt);
    for (unsigned int i = 0; i < ambient_dimension; i++) {
      lower[i] = (p[i] < lower[i]) ? p[i] : lower[i];
      upper[i] = (p[i] > upper[i]) ? p[i] : upper[i];
    }
  }

  double domain_size = 1;

  for (std::size_t i = 0; i < ambient_dimension; i++)
    domain_size *= CGAL::to_double(upper[i] - lower[i]);

  // named parameter for r_max
  double r_max = CGAL::to_double(parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_radius), 2 * internal::get_maximum_radius(dimension, number_of_points, domain_size)));
  double r_min = CGAL::to_double(weight_limiting ? internal::get_minimum_radius(points.size(), number_of_points, beta, gamma, r_max) : 0);

  auto weight_functor = parameters::choose_parameter(parameters::get_parameter(np, internal_np::weight_function), internal::Weight_functor<Point>(r_min, alpha));

  std::size_t heap_size = points.size();
  std::vector<Point> tiling_points;

  IPM ipm(points, tiling_points, point_map);

  auto tile_point = [&tiling_points, &lower, &upper, &r_max](const Point& p, std::size_t dim = 0) {
    auto do_tiling = [&tiling_points, &lower, &upper, &r_max](const auto& self, const Point& p, std::size_t dim) -> void {
      auto it = p.cartesian_begin();

      if (lower[int(dim)] > (*(it + dim) - r_max)) {
        FT v = 2 * lower[int(dim)] - (*(it + dim));
        Point p2;
        internal::copy_and_replace(p, p2, dim, v);
        tiling_points.emplace_back(p2);
        if (dim + 1 < CGAL::Ambient_dimension<Point>::value)
          self(self, tiling_points.back(), dim + 1);
      }

      if (upper[int(dim)] < (*(it + dim) + r_max)) {
        FT v = 2 * upper[int(dim)]- (*(it + dim));
        Point p2;
        internal::copy_and_replace(p, p2, dim, v);
        tiling_points.emplace_back(p2);
        if (dim + 1 < CGAL::Ambient_dimension<Point>::value)
          self(self, tiling_points.back(), dim + 1);
      }

      if (dim + 1 < CGAL::Ambient_dimension<Point>::value)
        self(self, p, dim + 1);
      };
    do_tiling(do_tiling, p, dim);
    };

  if (tiling) {
    for (std::size_t i = 0;i<points.size();i++)
      tile_point(get(point_map, *(points.begin() + i)));
  }

  // Tiling requires either to copy the point range or to use a special Index_property_map with a second array of points to keep the input range const.
  Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size() + tiling_points.size()), Splitter(), Search_traits(ipm));
  tree.build();

  std::vector<std::size_t> heap(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size()));
  std::vector<std::size_t> heap_pos(points.size());

  std::size_t target_points = number_of_points;

  do {
    // Computing weights
    std::vector<double> weights(points.size(), 0);
    std::vector<std::size_t> res;
    for (std::size_t i = 0; i < heap_size; i++) {
      const Point& p = get(ipm, heap[i]);
      res.clear();
      tree.search(std::back_inserter(res), Fuzzy_sphere(p, r_max, 0, Search_traits(ipm)));

      for (std::size_t n = 0; n < res.size(); n++) {
        if (i == res[n] || res[n] >= heap_pos.size() || heap_pos[res[n]] >= heap_size)
          continue;

        const Point &p2 = get(ipm, res[n]);
        double d2 = CGAL::to_double(squared_distance(p, p2));
        weights[i] += weight_functor(p, p2, d2, r_max);
      }
    }

    auto weight_order = [&](std::size_t a, std::size_t b) {
      return weights[a] < weights[b];
      };

    std::make_heap(heap.begin(), heap.begin() + heap_size, weight_order);

    // inverse mapping
    for (std::size_t i = 0; i < heap.size(); i++)
      heap_pos[heap[i]] = i;

    while (heap_size > target_points) {
      std::size_t i = heap[0];

      internal::pop_heap(heap, heap_pos, heap_size, weights);

      // Adjust weights of neighbors
      res.clear();
      const Point& p = get(point_map, *(points.begin() + i));
      tree.search(std::back_inserter(res), Fuzzy_sphere(p, r_max, 0, Search_traits(ipm)));

      for (std::size_t n = 0; n < res.size(); n++) {
        if (i == res[n] || res[n] >= points.size() || heap_pos[res[n]] >= heap_size)
          continue;

        const Point &p2 = get(point_map, *(points.begin() + res[n]));
        double d2 = CGAL::to_double(squared_distance(p, p2));

        weights[res[n]] -= weight_functor(p2, p, d2, r_max);

        internal::move_down(heap, heap_pos, heap_size, weights, heap_pos[res[n]]);
      }
    }

    CGAL_assertion(heap_size == target_points);

    if (progressive) {
      target_points = target_points>>1;
      r_max = r_max * std::pow(2.0, 1.0 / 3.0);
    }
  } while (progressive && target_points >= 3);

  for (std::size_t i = 0; i < number_of_points; i++)
    output++ = get(ipm, heap[i]);
}

} // namespace CGAL

#endif // CGAL_POISSON_ELIMINATE_H
