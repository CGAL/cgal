#ifndef CGAL_SURFACE_MESH_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
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


#define CGAL_SURFACE_MESH_SEGMENTATION_EXPECTATION_MAXIMIZATION_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>
#include <CGAL/assertions.h>
#include <CGAL/Random.h>

#define CGAL_DEFAULT_MAXIMUM_ITERATION 10u
#define CGAL_DEFAULT_NUMBER_OF_RUN 15u
#define CGAL_DEFAULT_THRESHOLD 1e-3

#define CGAL_DEFAULT_SEED 1340818006

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{

/**
 * @brief Expectation maximization algorithm for GMM fitting.
 * @see Gaussian_center
 */
class Expectation_maximization
{
private:
  /**
   * @brief Represents centers in Expectation Maximization algorithm.
   * @see Expectation_maximization
   */
  class Gaussian_center
  {
  public:
    double mean;
    double deviation;
    double mixing_coefficient;

    Gaussian_center(): mean(0), deviation(0), mixing_coefficient(1.0) {
    }
    Gaussian_center(double mean, double deviation = 0.0,
                    double mixing_coefficient = 0.0)
      : mean(mean), deviation(deviation), mixing_coefficient(mixing_coefficient) {
    }
    operator double() const {
      return mean;
    }
    /**
     * Probability density function (pdf).
     * Note that result is not devided to \f$ \sqrt {2\pi}  \f$ , since it does not effect EM algorithm.
     * @param x data
     * @return pdf result (without dividing \f$ \sqrt {2\pi}  \f$)
     */
    double probability(double x) const {
      if(deviation == 0.0) {
        return x == mean ? 1.0 : 0.0;
      }

      double e_over = -0.5 * std::pow((x - mean) / deviation, 2);
      return exp(e_over) / deviation;
    }
    /**
     * Multiplies pdf result and mixing coefficient of the center.
     * @param x data
     * @return result of the multiplication.
     */
    double probability_with_coef(double x) const {
      return probability(x) * mixing_coefficient;
    }
    /** A comparator for sorting centers in ascending order. */
    bool operator < (const Gaussian_center& center) const {
      return mean < center.mean;
    }
  };
public:
  /** Options for initial center placement. */
  enum Initialization_types {
    RANDOM_INITIALIZATION, /**< place initial centers randomly */
    PLUS_INITIALIZATION,   /**< place initial centers using k-means++ algorithm */
    K_MEANS_INITIALIZATION /**< run k-means clustering and use result of it as initial center positions */
  };

  double final_likelihood;
private:
  std::vector<Gaussian_center>      centers;
  std::vector<double>               points;
  std::vector<std::vector<double> > responsibility_matrix;

  double threshold;
  std::size_t    maximum_iteration;

  Initialization_types init_type;

  CGAL::Random random;
public:
  /**
   * @pre @a number_of_centers should be positive
   * @pre size of @a data should be no smaller than number_of_centers
   *
   * Constructs structures and runs the algorithm.
   *
   * If @a init_type is either RANDOM_INITIALIZATION or PLUS_INITIALIZATION,
   * then EM algorithm is repeated @a number_of_runs times, and the result which has maximum likelihood is kept.
   * Otherwise (i.e. init_type is K_MEANS_INITIALIZATION) EM algorithm is just run one time.
   * @param number_of_centers
   * @param data
   * @param init_type option for initial center placement.
   * @param number_of_runs number of times to repeat EM algorithm
   * @param threshold minimum allowed improvement on likelihood between iterations
   * @param maximum_iteration maximum allowed iteration in a single EM algorithm call
   */
  Expectation_maximization(std::size_t number_of_centers,
                           const std::vector<double>& data,
                           Initialization_types init_type = PLUS_INITIALIZATION,
                           std::size_t number_of_runs = CGAL_DEFAULT_NUMBER_OF_RUN,
                           double threshold = CGAL_DEFAULT_THRESHOLD,
                           std::size_t maximum_iteration = CGAL_DEFAULT_MAXIMUM_ITERATION )
    :
    final_likelihood(-(std::numeric_limits<double>::max)()), points(data),
    responsibility_matrix(std::vector<std::vector<double> >(number_of_centers,
                          std::vector<double>(points.size()))),
    threshold(threshold),
    maximum_iteration(maximum_iteration),
    init_type(init_type),
    random(CGAL_DEFAULT_SEED) {
    CGAL_assertion(data.size() >= number_of_centers
                   && "Number of centers can not be more than number of data.");

    // For initialization with k-means, with one run
    if(init_type == K_MEANS_INITIALIZATION) {
      K_means_clustering k_means(number_of_centers, data,
                                 K_means_clustering::PLUS_INITIALIZATION,
                                 number_of_runs, maximum_iteration);

      k_means.fill_with_centers(centers);
      calculate_initial_mixing_and_deviation();

      calculate_clustering();
    }
    // For initialization with random center selection, with multiple run
    else {
      calculate_clustering_with_multiple_run(number_of_centers, number_of_runs);
    }
    sort(centers.begin(), centers.end());
  }

  /**
   * Fills data_center by the id of the center which has maximum responsibility.
   * @param[out] data_centers
   */
  void fill_with_center_ids(std::vector<std::size_t>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<double>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      double max_likelihood = 0.0;
      std::size_t max_center = (std::numeric_limits<std::size_t>::max)(),
                  center_counter = 0;
      for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it, ++center_counter) {
        double likelihood = center_it->probability_with_coef(*point_it);
        if(max_likelihood < likelihood) {
          max_likelihood = likelihood;
          max_center = center_counter;
        }
      }
      CGAL_assertion( max_center!=(std::numeric_limits<std::size_t>::max)() );
      data_centers.push_back(max_center);
    }
  }

  /**
   * Fills probabilities[center][point] by responsibility of the center on the point.
   * @param[out] probabilities
   */
  void fill_with_probabilities(std::vector<std::vector<double> >& probabilities) {
    probabilities = std::vector<std::vector<double> >
                    (centers.size(), std::vector<double>(points.size()));
    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      double total_probability = 0.0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double probability = centers[center_i].probability_with_coef(points[point_i]);
        total_probability += probability;
        probabilities[center_i][point_i] = probability;
      }
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        probabilities[center_i][point_i] /= total_probability;
      }
    }
  }

private:
  /**
   * Calculates deviation for each center.
   * Initial deviation of a center is equal to deviation of the points whose closest center is the current center.
   */
  void calculate_initial_mixing_and_deviation() {
    // assign same mixing coef for each cluster
    for(std::vector<Gaussian_center>::iterator it = centers.begin();
        it != centers.end(); ++it) {
      it->mixing_coefficient = 1.0 / centers.size();
    }

    // calculate deviation
    std::vector<std::size_t> member_count(centers.size(), 0);
    for(std::vector<double>::iterator it = points.begin(); it!= points.end();
        ++it) {
      std::size_t closest_center = 0;
      double min_distance = std::abs(centers[0].mean - *it);
      for(std::size_t i = 1; i < centers.size(); ++i) {
        double distance = std::abs(centers[i].mean - *it);
        if(distance < min_distance) {
          min_distance = distance;
          closest_center = i;
        }
      }
      member_count[closest_center]++;
      centers[closest_center].deviation += min_distance * min_distance;
    }
    for(std::size_t i = 0; i < centers.size(); ++i) {
      if(member_count[i] == 0) {
        CGAL_assertion( false &&
                        "There is a cluster which does not contain any points, it will not cause an error but associated probabilites to this cluster will be 0.");
      } else {
        centers[i].deviation = std::sqrt(centers[i].deviation / member_count[i]);
      }
    }
  }

  /**
   * Initializes centers by choosing random points from data.
   * @param number_of_centers
   */
  void initiate_centers_randomly(std::size_t number_of_centers) {
    centers.clear();
    Selector().forgy_initialization(number_of_centers, points, centers, random);

    calculate_initial_mixing_and_deviation();
  }

  /**
   * Initializes centers by using K-means++ algorithm.
   * Probability of a point to become a center is proportional to its squared distance to the closest center.
   * @param number_of_centers
   */
  void initiate_centers_plus_plus(std::size_t number_of_centers) {
    centers.clear();
    Selector().plus_plus_initialization(number_of_centers, points, centers, random);

    calculate_initial_mixing_and_deviation();
  }

  //Main steps of EM algorithm

  /**
   * Corresponds to M step.
   * Recalculates parameters of the centers using current responsibility matrix.
   */
  void calculate_parameters() {
    for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
      // Calculate new mean
      double new_mean = 0.0, total_membership = 0.0;
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double membership = responsibility_matrix[center_i][point_i];
        new_mean += membership * points[point_i];
        total_membership += membership;
      }
      new_mean /= total_membership;

      // Calculate new deviation
      double new_deviation = 0.0;
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double membership = responsibility_matrix[center_i][point_i];
        new_deviation += membership * std::pow(points[point_i] - new_mean, 2);
      }
      new_deviation = std::sqrt(new_deviation/total_membership);

      // Assign new parameters
      centers[center_i].mixing_coefficient = total_membership / points.size();
      centers[center_i].deviation = new_deviation;
      centers[center_i].mean = new_mean;
    }
  }

  /**
   * Corresponds to both E step and likelihood step.
   * Calculates log-likelihood, and responsibility matrix using current center parameters.
   * @return log-likelihood
   */
  double calculate_likelihood() {
    // The trick (merely a trick) is while calculating log-likelihood, we also refresh responsibility matrix,
    // so that in next iteration we do not have to calculate matrix again.

    double likelihood = 0.0;
    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      double total_membership = 0.0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double membership = centers[center_i].probability_with_coef(points[point_i]);
        total_membership += membership;
        responsibility_matrix[center_i][point_i] = membership;
      }
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        responsibility_matrix[center_i][point_i] /= total_membership;
      }
      likelihood += log(total_membership);
    }
    return likelihood;
  }

  /**
   * One iteration of EM algorithm. Includes E-step, M-step and likelihood calculation.
   * @param first_iteration
   * @return log-likelihood
   * @see calculate_likelihood() for E-step and likelihood calculation, calculate_parameters() for M-step
   */
  double iterate(bool first_iteration) {
    // E-step
    // we call calculate_likelihood for E-step in first iteration because
    // at first iteration, E-step is not done since calculate_likelihood() is not called yet.
    if(first_iteration) {
      calculate_likelihood();
    }

    // M-step
    calculate_parameters();

    // Likelihood step and also E-step for next iteration
    return calculate_likelihood(); // calculates likelihood and -also- refreshes responsibility matrix,
    // so that we do not have to calculate it in next iteration.
  }

  /**
   * Main entry point for EM algorithm.
   * Iterates until convergence occurs (i.e. likelihood - prev_likelihood < threshold * std::abs(likelihood))
   * or maximum iteration limit is reached.
   * @see iterate()
   */
  double calculate_clustering() {
    double likelihood = -(std::numeric_limits<double>::max)(), prev_likelihood;
    std::size_t iteration_count = 0;
    double is_converged = false;
    while(!is_converged && iteration_count++ < maximum_iteration) {
      prev_likelihood = likelihood;
      likelihood = iterate(iteration_count == 1);
      double progress = likelihood - prev_likelihood;
      is_converged = progress < threshold * std::abs(likelihood);
    }
    if(final_likelihood < likelihood) {
      final_likelihood = likelihood;
    }
    return likelihood;
  }

  /**
   * Calls calculate_clustering() @a number_of_run times,
   * and keeps the result which has maximum likelihood.
   * @param number_of_centers
   * @param number_of_run
   * @see calculate_clustering()
   */
  void calculate_clustering_with_multiple_run(std::size_t number_of_centers,
      std::size_t number_of_run) {
    std::vector<Gaussian_center> max_centers;

    while(number_of_run-- > 0) {
      init_type == RANDOM_INITIALIZATION ? initiate_centers_randomly(
        number_of_centers)
      : initiate_centers_plus_plus(number_of_centers);

      double likelihood = calculate_clustering();
      if(likelihood == final_likelihood) {
        max_centers = centers;
      }
    }
    centers = max_centers;
  }
};
}//namespace internal
/// @endcond
}//namespace CGAL
#undef CGAL_DEFAULT_SEED
#undef CGAL_DEFAULT_MAXIMUM_ITERATION
#undef CGAL_DEFAULT_THRESHOLD
#undef CGAL_DEFAULT_NUMBER_OF_RUN
#endif //CGAL_SURFACE_MESH_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
