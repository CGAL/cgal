#ifndef CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Ilker O. Yaz


#define CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H

#include <CGAL/license/Surface_mesh_segmentation.h>


#include <vector>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>

#include <CGAL/assertions.h>
#include <CGAL/Random.h>

#define CGAL_DEFAULT_MAXIMUM_ITERATION 10u
#define CGAL_DEFAULT_NUMBER_OF_RUN 15u
#define CGAL_DEFAULT_SEED 1340818006

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{

/**
 * Class providing initialization functionality: Forgy initialization, k++ initialization
 *
 * TODO: rand() might be changed with a generator for larger 'range'
 */
class Selector
{
public:
  /**
   * Selects random samples from `points` and put them into `centers`. Does not select one sample more than once as center.
   * T2 should be constructable by T1
   *
   * Implementation note: it is a variant of Floyd generator, and has uniform distribution
   * where k = number of centers = complexity is O(k log k), and mem overhead is O(k)
   *
   * I also left previous implementation below, it might be useful where number of centers close to number of points
   */
  template<class T1, class T2>
  void forgy_initialization(std::size_t number_of_centers,
                            const std::vector<T1>& points, std::vector<T2>& centers, CGAL::Random& random) {
    centers.reserve(number_of_centers);
    std::set<std::size_t> selected;

    for(std::size_t i = 0; i < number_of_centers; ++i) {
      std::size_t random_range = points.size() - number_of_centers +
                                 i; // activate one more element in each iteration for as selectable
      std::size_t random_index = random(random_range + 1); // [0, random_range];

      std::pair<std::set<std::size_t>::iterator, bool> random_index_unique =
        selected.insert(random_index);
      if(!random_index_unique.second) {
        random_index_unique = selected.insert(
                                random_range); // random_range can not be inside random_index_unique,
        CGAL_assertion(
          random_index_unique.second);          // since we just increment it by one
      }

      centers.push_back(points[ *(random_index_unique.first) ]);
    }
  }

  // To future reference, I also left prev implementation which is a variant of Fisher–Yates shuffle, however to keep `points` intact I use another vector to
  // store and swap indices.
  // where n = number of points; complexity = O(n), memory overhead = O(n)
  /*
  template<class T1, class T2>
  void forgy_initialization(std::size_t number_of_centers, const std::vector<T1>& points, std::vector<T2>& centers)
  {
      std::vector<std::size_t> indices(points.size()); // it is required to not to swap points
      for(std::size_t i = 0; i < points.size(); ++i) { indices[i] = i; }

      centers.reserve(number_of_centers);
      for(std::size_t i = 0; i < number_of_centers; ++i)
      {
          std::size_t random_index = rand() % (points.size() - i); // select a random index between 0 and not selected sample count
          centers.push_back(points[ indices[random_index] ]);

          std::swap(indices[random_index], indices[points.size() -i -1]); // not to select again, swap at the current-end
      }
  }
  */

  /**
   * Selects samples from `points` using K-means++ algorithm and put them into `centers`.
   * Does not select one sample more than once as center.
   * Probability of a point to become a center is proportional to its squared distance to the closest center.
   *
   * T2 should be constructable by T1, and both T1 and T2 should have conversion operator to double.
   */
  template<class T1, class T2>
  void plus_plus_initialization(std::size_t number_of_centers,
                                const std::vector<T1>& points, std::vector<T2>& centers, CGAL::Random& random) {
    centers.reserve(number_of_centers);

    std::vector<double> distance_square(points.size(),
                                        (std::numeric_limits<double>::max)());
    std::vector<double> distance_square_cumulative(points.size());

    // distance_square stores squared distance to the closest center for each point.
    // say, "distance_square" ->            [ 0.1, 0.2, 0.3, 0.0, 0.2 ... ]
    // then distance_square_cumulative ->   [ 0.1, 0.3, 0.6, 0.6, 0.8 ... ]
    std::size_t initial_index = random(points.size()); // [0, points size)
    centers.push_back(points[initial_index]);

    for(std::size_t i = 1; i < number_of_centers; ++i) {
      double cumulative_distance_square = 0.0;
      // distance_square holds closest distance that points have, so just test new coming center (i.e. centers.back())
      for(std::size_t j = 0; j < points.size(); ++j) {
        double new_distance = std::pow(centers.back() - points[j], 2);
        if(new_distance < distance_square[j]) {
          distance_square[j] = new_distance;
        }
        cumulative_distance_square += distance_square[j];
        distance_square_cumulative[j] = cumulative_distance_square;
      }

      // check whether there is a way to select a new center
      double total_probability = distance_square_cumulative.back();
      CGAL_assertion( total_probability != 0.0 &&
                      "Squared distances of all elements to their closest centers are 0.0, there is no way to select a new center.");

      double random_ds = random.get_double(0.0,
                                           total_probability); // [0.0, total_probability)

      // this can not select end(), since random_ds < total_probability (i.e. distance_square_cumulative.back())
      // this can not select an already selected item since either (by considering that upper bounds returns greater)
      //  - aready selected item is at 0, and its value is 0.0
      //  - or its value is equal to value of previous element
      std::size_t selection_index = std::upper_bound(
                                      distance_square_cumulative.begin(), distance_square_cumulative.end(), random_ds)
                                    - distance_square_cumulative.begin();

      centers.push_back(points[selection_index]);
    }
  }
};

class K_means_center;

/**
  * @brief Represents points in k-means algorithm.
  * @see K_means_center, K_means_clustering
  */
class K_means_point
{
public:
  double data;      /**< Location of the point */
  std::size_t    center_id; /**< Closest center to the point */
  K_means_point(double data,
                std::size_t center_id = (std::numeric_limits<std::size_t>::max)())
    : data(data), center_id(center_id) {
  }

  operator double() const {
    return data;
  }

  bool calculate_new_center(std::vector<K_means_center>& centers);
};

/**
  * @brief Represents centers in k-means algorithm.
  * @see K_means_point, K_means_clustering
  */
class K_means_center
{
public:
  double mean; /**< Mean of the center */
private:
  double new_mean;
  std::size_t    new_number_of_points;
  bool   empty;

public:

  K_means_center(double mean): mean(mean), new_mean(0.0), new_number_of_points(0),
    empty(true) {
  }

  K_means_center(const K_means_point& point) : mean(point.data), new_mean(0.0),
    new_number_of_points(0), empty(true) {
  }

  operator double() const {
    return mean;
  }

  /**
    * Called by a point for adding itself to mean of the center (i.e. registering).
    * @param data location of the point
    */
  void add_point(double data) {
    ++new_number_of_points;
    new_mean += data;
  }
  /**
    * Called after every point is registered to its closest center
    * @see add_point()
    */
  void calculate_mean() {
    if(new_number_of_points ==
        0) { // just return and this will leave `mean` as it is,
      empty = true;                 // `new_mean` is 0.0 since `new_number_of_points` is 0
      return;
    }

    empty = false;
    mean = new_mean / new_number_of_points;
    new_number_of_points = 0;
    new_mean = 0.0;
  }

  bool is_empty() {
    return empty;
  }

  /** A comparator for sorting centers in ascending order. */
  bool operator < (const K_means_center& center) const {
    return mean < center.mean;
  }
};

/**
* Finds closest center and adds itself to the closest center's mean.
* @param centers available centers
* @return true if #center_id is changed (i.e. new center != previous center)
*/
inline bool K_means_point::calculate_new_center(std::vector<K_means_center>&
    centers)
{
  std::size_t new_center_id = 0;
  double min_distance = std::abs(centers[0].mean - data);
  for(std::size_t i = 1; i < centers.size(); ++i) {
    double new_distance = std::abs(centers[i].mean - data);
    if(new_distance < min_distance) {
      new_center_id = i;
      min_distance = new_distance;
    }
  }
  bool is_center_changed = (new_center_id != center_id);
  center_id = new_center_id;

  centers[center_id].add_point(data);
  return is_center_changed;
}


/**
 * @brief K-means clustering algorithm.
 * @see K_means_point, K_means_center
 */
class K_means_clustering
{
public:
  /** Types of algorithms for random center selection. */
  enum Initialization_types {
    RANDOM_INITIALIZATION, /**< place initial centers randomly */
    PLUS_INITIALIZATION    /**< place initial centers using k-means++ algorithm */
  };

private:
  std::vector<K_means_center> centers;
  std::vector<K_means_point>  points;
  std::size_t  maximum_iteration;

  Initialization_types init_type;

  CGAL::Random random;
public:
  /**
   * @pre @a number_of_centers should be positive
   * @pre size of @a data should be no smaller than number_of_centers
   *
   * Constructs structures and runs the algorithm.
   * K-means algorithm is repeated number_of_run times, and the result which has minimum within cluster error is kept.
   * @param number_of_centers
   * @param data
   * @param init_type initialization type for random center selection
   * @param number_of_run number of times to repeat k-means algorithm
   * @param maximum_iteration maximum allowed iteration in a single k-means algorithm call
   *
   *
   */
  K_means_clustering(std::size_t number_of_centers,
                     const std::vector<double>& data,
                     Initialization_types init_type = PLUS_INITIALIZATION,
                     std::size_t number_of_run = CGAL_DEFAULT_NUMBER_OF_RUN,
                     std::size_t maximum_iteration = CGAL_DEFAULT_MAXIMUM_ITERATION)
    :
    points(data.begin(), data.end()),
    maximum_iteration(maximum_iteration),
    init_type(init_type),
    random(CGAL_DEFAULT_SEED) {
    CGAL_precondition(data.size() >= number_of_centers
                      && "Number of centers can not be more than number of data.");

    calculate_clustering_with_multiple_run(number_of_centers, number_of_run);
    sort(centers.begin(), centers.end());
  }

  /**
   * Fills data_center by the id of the closest center for each point.
   * @param[out] data_centers
   */
  void fill_with_center_ids(std::vector<std::size_t>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      point_it->calculate_new_center(centers); // find closest center
      data_centers.push_back(point_it->center_id);
    }
  }

  template<class T>
  void fill_with_centers(std::vector<T>& center_means) {
    for(std::vector<K_means_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      center_means.push_back(center_it->mean);
    }
  }

private:
  /**
   * Initializes centers by choosing random points from data. Does not select one points more than once as center.
   * @param number_of_centers
   */
  void initiate_centers_randomly(std::size_t number_of_centers) {
    centers.clear();
    Selector().forgy_initialization(number_of_centers, points, centers, random);
  }

  /**
   * Initializes centers by using K-means++ algorithm.
   * Probability of a point to become a center is proportional to its squared distance to the closest center.
   * @param number_of_centers
   */
  void initiate_centers_plus_plus(std::size_t number_of_centers) {
    centers.clear();
    Selector().plus_plus_initialization(number_of_centers, points, centers, random);
  }

  /**
   * One iteration of k-means algorithm.
   * @return true if any closest center to a point is changed
   */
  bool iterate() {
    bool any_center_changed = false;
    // For each point, calculate its new center
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      bool center_changed = point_it->calculate_new_center(centers);
      any_center_changed |= center_changed;
    }
    // For each center, calculate its new mean
    for(std::vector<K_means_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      center_it->calculate_mean();
    }
    return any_center_changed;
  }

  /**
   * Main entry point for k-means algorithm.
   * Iterates until convergence occurs (i.e. no point changes its center) or maximum iteration limit is reached.
   */
  void calculate_clustering() {
    std::size_t iteration_count = 0;
    bool any_center_changed = true;
    while(any_center_changed && iteration_count++ < maximum_iteration) {
      any_center_changed = iterate();
    }
  }

  /**
   * Calls calculate_clustering() @a number_of_run times,
   * and keeps the result which has minimum within cluster error.
   * @param number_of_centers
   * @param number_of_run
   * @see calculate_clustering(), within_cluster_sum_of_squares()
   */
  void calculate_clustering_with_multiple_run(std::size_t number_of_centers,
      std::size_t number_of_run) {
    std::vector<K_means_center> min_centers;
    double error = (std::numeric_limits<double>::max)();
    while(number_of_run-- > 0) {
      init_type == RANDOM_INITIALIZATION ? initiate_centers_randomly(
        number_of_centers)
      : initiate_centers_plus_plus(number_of_centers);
      calculate_clustering();
      double new_error = within_cluster_sum_of_squares();
      if(error > new_error) {
        error = new_error;
        min_centers = centers;
      }
    }
    // Note that current center-ids in points are not valid.
    // But they are recalculated when asked (in fill_with_center_ids())
    centers = min_centers;
  }

  /**
   * Sum of squared distances between each point and the closest center to it.
   */
  double within_cluster_sum_of_squares() const {
    double sum = 0.0;
    for(std::vector<K_means_point>::const_iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      sum += std::pow(centers[point_it->center_id].mean - point_it->data, 2);
    }
    return sum;
  }
};
}//namespace internal
/// @endcond
}//namespace CGAL
#undef CGAL_DEFAULT_SEED
#undef CGAL_DEFAULT_MAXIMUM_ITERATION
#undef CGAL_DEFAULT_NUMBER_OF_RUN

#endif //CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H
