// Copyright (c) 2024 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H
#define CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#ifdef CGAL_USE_CY
#include <cyVector.h>
#include <cySampleElim.h>
#endif


#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {
template <class FT>
FT get_maximum_radius(std::size_t dimensions, std::size_t sample_size, FT domain_size) {
  FT sampleArea = domain_size / (FT)sample_size;
  FT r_max;
  switch (dimensions) {
  case 2: r_max = CGAL::sqrt(sampleArea / (FT(2) * CGAL::sqrt(FT(3)))); break;
  case 3: r_max = std::pow(sampleArea / (FT(4) * CGAL::sqrt(FT(2))), FT(1) / FT(3)); break;
  default:
    FT c;
    std::size_t d_start;
    if ((dimensions & 1)) { c = FT(2);      d_start = 3; }
    else { c = CGAL_PI; d_start = 4; }
    for (std::size_t d = d_start; d <= dimensions; d += 2) c *= FT(2) * CGAL_PI / FT(d);
    r_max = std::pow(sampleArea / c, FT(1) / FT(dimensions));
    break;
  }
  return r_max;
}

template <class FT>
FT get_minimum_radius(std::size_t inputSize, std::size_t outputSize, FT beta, FT gamma, FT r_max) {
  FT ratio = FT(outputSize) / FT(inputSize);
  return r_max * (1 - std::pow(ratio, gamma)) * beta;
}

template<class Point, class FT = CGAL::Kernel_traits<Point>::FT, unsigned int dim = CGAL::Ambient_dimension<Point>::value>
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

  out = construct<Point>(tmp);
}

template<class PointRange, class PointMap>
class Indexed_point_map {
public:
  using value_type = typename boost::property_traits<PointMap>::value_type;
  using reference = const value_type&;
  using key_type = std::size_t;
  using category = boost::readable_property_map_tag;

private:
  const PointRange& range;
  const std::vector<value_type>& tiling_points;
  const PointMap& point_map;

public:
  Indexed_point_map(const PointRange& input, const std::vector<value_type>& tiling, const PointMap& point_map)
    : range(input), tiling_points(tiling), point_map(point_map) {}

  reference operator[](key_type k) const {
    CGAL_assertion(k < (range.size() + tiling_points.size()));
    if (range.size() > k)
      return get(point_map, range[k]);
    else
      return tiling_points[k - range.size()];
  }

  friend reference get(const Indexed_point_map& ppmap, key_type k) {
    CGAL_assertion(k < (ppmap.range.size() + ppmap.tiling_points.size()));
    if ((k < ppmap.range.size()))
      return get(ppmap.point_map, ppmap.range[k]);
    else
      return ppmap.tiling_points[k - ppmap.range.size()];
  }
};


template<class Point_3>
class Weight_functor {
public:
  using FT = typename Kernel_traits<Point_3>::Kernel::FT;
  Weight_functor(FT r_min = 0, FT alpha = 8) : r_min(r_min), alpha(alpha) {}

  FT operator()(const Point_3&, const Point_3&, FT d2, FT r_max) {
    FT d = CGAL::sqrt(d2);
    if (d < r_min) d = r_min;
    return std::pow(FT(1) - d / r_max, alpha);
  }

private:
  FT r_min;
  FT alpha;
};

template<class FT>
void move_down(std::vector<std::size_t>& heap, std::vector<std::size_t>& heap_pos, std::size_t heap_size, std::vector<FT>& weights, std::size_t idx) {
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

template<class FT>
void pop_heap(std::vector<std::size_t>& heap, std::vector<std::size_t>& heap_pos, std::size_t &heap_size, std::vector<FT>& weights) {
  std::swap(heap.front(), heap[heap_size - 1]);
  heap_pos[heap.front()] = 0;
  heap_pos[heap[heap_size - 1]] = heap_size - 1;

  heap_size--;

  move_down(heap, heap_pos, heap_size, weights, 0);
}

/*
template<class FT>
void check_heap(std::vector<std::size_t>& heap, std::size_t heap_size, std::vector<FT>& weights) {
  std::size_t idx = 0;
  std::size_t child = 2 * idx + 1;

  while (child + 1 < heap_size) {
    if (weights[heap[idx]] < weights[heap[child]])
      std::cout << std::endl;
    if (weights[heap[idx]] < weights[heap[child + 1]])
      std::cout << std::endl;

    idx++;
    child = 2 * idx + 1;
  }

  if (child < heap_size) {
    if (weights[heap[idx]] < weights[heap[child]])
      std::cout << std::endl;
  }
}
*/

}

template<class PointRange, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
void poisson_eliminate(PointRange points, std::size_t number_of_points, OutputIterator out, const NamedParameters& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters>;
  using PointMap = typename NP_helper::Point_map;
  using Point = typename boost::property_traits<PointMap>::value_type;
  using GeomTraits = typename Kernel_traits<Point>::Kernel;
  using FT = typename GeomTraits::FT;
  using IPM = internal::Indexed_point_map<PointRange, PointMap>;
  PointMap point_map = NP_helper::get_point_map(points, np);

  CGAL::Real_timer timer, heap_t;

  using Search_traits = CGAL::Search_traits_adapter<std::size_t, IPM, CGAL::Search_traits_3<GeomTraits>>;
  //using Splitter = CGAL::Sliding_midpoint<Search_traits>;
  using Splitter = CGAL::Sliding_midpoint<Search_traits>;
  using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;
  using Tree = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true, CGAL::Tag_true>;

  // Needs a multi-dimension solution
  Bbox_3 bb = bbox_3(make_transform_iterator_from_property_map(points.begin(), point_map),
    make_transform_iterator_from_property_map(points.end(), point_map));

  FT domain_size = (bb.xmax() - bb.xmin()) * (bb.ymax() - bb.ymin()) * (bb.zmax() - bb.zmin());

  // named parameters for alpha, beta and gamma
  const FT alpha = 8;
  const FT beta = 0.65;
  const FT gamma = 1.5;
  // named parameters for weight limiting
  const bool weight_limiting = parameters::choose_parameter(parameters::get_parameter(np, internal_np::weight_limiting), true);
  // named parameter for progressive
  const bool progressive = parameters::choose_parameter(parameters::get_parameter(np, internal_np::progressive), false);
  // named parameter for tiling
  const bool tiling = parameters::choose_parameter(parameters::get_parameter(np, internal_np::tiling), false);
  const unsigned int ambient_dimension = CGAL::Ambient_dimension<Point>::value;
  const unsigned int dimension = parameters::choose_parameter(parameters::get_parameter(np, internal_np::dimension), ambient_dimension);

  // named parameter for r_max
  FT r_max = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_radius), 2 * internal::get_maximum_radius(dimension, number_of_points, domain_size));
  FT r_min = (weight_limiting ? internal::get_minimum_radius(points.size(), number_of_points, beta, gamma, r_max) : 0);

  auto weight_functor = parameters::choose_parameter(parameters::get_parameter(np, internal_np::weight_functor), internal::Weight_functor<Point>(r_min, alpha));

  std::size_t heap_size = points.size();
  std::size_t u = 0;
  std::vector<Point> tiling_points;

  IPM ipm(points, tiling_points, point_map);

  auto tile_point = [&tiling_points, &bb, &r_max, &ambient_dimension](const Point& p, std::size_t dim = 0) {
    auto do_tiling = [&tiling_points, &bb, &r_max, &ambient_dimension](const auto& self, const Point& p, std::size_t dim) -> void {
      auto it = p.cartesian_begin();

      if (bb.min_coord(int(dim)) > (*(it + dim) - r_max)) {
        FT v = 2 * bb.min_coord(int(dim)) - (*(it + dim));
        Point p2;
        internal::copy_and_replace(p, p2, dim, v);
        tiling_points.emplace_back(p2);
        if (dim + 1 < ambient_dimension)
          self(self, tiling_points.back(), dim + 1);
      }

      if (bb.max_coord(int(dim)) < (*(it + dim) + r_max)) {
        FT v = 2 * bb.max_coord(int(dim)) - (*(it + dim));
        Point p2;
        internal::copy_and_replace(p, p2, dim, v);
        tiling_points.emplace_back(p2);
        if (dim + 1 < ambient_dimension)
          self(self, tiling_points.back(), dim + 1);
      }

      if (dim + 1 < ambient_dimension)
        self(self, p, dim + 1);
      };
    do_tiling(do_tiling, p, dim);
    };

  if (tiling) {
    for (std::size_t i = 0;i<points.size();i++)
      tile_point(get(point_map, points[i]));
  }

  // Tiling requires either to copy the point range or to use a special Index_property_map with a second array of points to keep the input range const.
  Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size() + tiling_points.size()), Splitter(), Search_traits(ipm));
  tree.build();

  std::vector<std::size_t> heap(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size()));
  std::vector<std::size_t> heap_pos(points.size());

  std::size_t target_points = number_of_points;

  do {
    // Computing weights
    std::vector<FT> weights(points.size(), 0);
    std::vector<std::size_t> res;
    for (std::size_t i = 0; i < heap_size; i++) {
      const Point& p = get(ipm, heap[i]);
      res.clear();
      tree.search(std::back_inserter(res), Fuzzy_sphere(p, r_max, 0, Search_traits(ipm)));

      for (std::size_t n = 0; n < res.size(); n++) {
        // Why heap_pos[res[n]] >= heap_size? it excludes points from former progressive steps?
        // How do I still consider tiling points?
        if (i == res[n] || res[n] >= heap_pos.size() || heap_pos[res[n]] >= heap_size)
          continue;

        const Point p2 = get(ipm, res[n]);
        FT d2 = (p - p2).squared_length();
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
      const Point& p = get(point_map, points[i]);
      tree.search(std::back_inserter(res), Fuzzy_sphere(p, r_max, 0, Search_traits(ipm)));

      for (std::size_t n = 0; n < res.size(); n++) {
        if (i == res[n] || res[n] >= points.size() || heap_pos[res[n]] >= heap_size)
          continue;

        const Point p2 = get(point_map, points[res[n]]);
        FT d2 = (p - p2).squared_length();

        weights[res[n]] -= weight_functor(p2, p, d2, r_max);

        internal::move_down(heap, heap_pos, heap_size, weights, heap_pos[res[n]]);
      }
    }

    CGAL_assertion(heap_size == target_points);

    if (progressive) {
      target_points = target_points>>1;
      r_max = r_max * FT(std::pow(2.0, (1.0) / double(3.0)));
    }
  } while (progressive && target_points >= 3);

  for (std::size_t i = 0; i < number_of_points; i++)
    out++ = get(ipm, heap[i]);
}
/*

template <class TriangleMesh, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
void poisson_eliminate(const TriangleMesh& sm, OutputIterator out, const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GeomTraits;
  typedef typename GeomTraits::Point_3                                            Point_3;

  Bbox_3 bb = bbox_3(sm.points().begin(), sm.points().end());
  cy::Vec3d bl(bb.xmin(), bb.ymin(), bb.zmin());
  cy::Vec3d tr(bb.xmax(), bb.ymax(), bb.zmax());

  std::vector<cy::Vec3d> inputPoints, outputPoints;
  std::vector<Point_3> points, output;

  // @todo write with a transform_iterator directly into inputPoints
  sample_triangle_mesh(sm,
                       std::back_inserter(points),
                       CGAL::parameters::number_of_points_on_faces(2* num_vertices(sm))
                         .do_sample_vertices(false)
                         .do_sample_edges(false));
  double area = CGAL::Polygon_mesh_processing::area(sm);

  for(int i = 0; i < points.size(); ++i){
    inputPoints.push_back(cy::Vec3d(to_double(points[i].x()), to_double(points[i].y()), to_double(points[i].z())));
  }

  outputPoints.resize(num_vertices(sm)/2);

  CGAL::IO::write_points("orig_poisson.xyz", points, CGAL::parameters::stream_precision(17));

  cy::WeightedSampleElimination< cy::Vec3d, double, 3, int > wse;
  wse.SetBoundsMin(bl);
  wse.SetBoundsMax(tr);
  bool isProgressive = true;

  double d_max = 2 * wse.GetMaxPoissonDiskRadius( 2, outputPoints.size(), area );

  poisson_eliminate2(points, outputPoints.size(), std::back_inserter(output), parameters::maximum_radius(d_max));
  CGAL::IO::write_points("my_poisson.xyz", output, CGAL::parameters::stream_precision(17));

  wse.Eliminate( inputPoints.data(), inputPoints.size(),
                 outputPoints.data(), outputPoints.size(),
                 isProgressive,
                 d_max, 2 );

  for (const cy::Vec3d& p : outputPoints){
    *out++ = Point_3(p.x, p.y, p.z);
  }
}*/

} // namespace Polygon_mesh_processing
} // namespace CGAL

//#endif // ifdef CGAL_USE_CY

#endif // CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H
