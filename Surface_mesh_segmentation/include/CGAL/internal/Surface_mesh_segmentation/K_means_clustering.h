#ifndef CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H
#define CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H

#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <algorithm>

#define CGAL_DEFAULT_MAXIMUM_ITERATION 10
#define CGAL_DEFAULT_NUMBER_OF_RUN 15
#define CGAL_DEFAULT_SEED 1340818006

namespace CGAL
{
namespace internal
{

/**
 * @brief K-means clustering algorithm.
 * @see K_means_point, K_means_center
 */
class K_means_clustering
{
// Nested classes
private:
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
    int    new_number_of_points;

  public:
    K_means_center(double mean): mean(mean), new_mean(0.0),
      new_number_of_points(0) {
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
      mean = new_mean / new_number_of_points;
      new_number_of_points = 0;
      new_mean = 0.0;
    }
    /** A comparator for sorting centers in ascending order. */
    bool operator < (const K_means_center& center) const {
      return mean < center.mean;
    }
  };

  /**
   * @brief Represents points in k-means algorithm.
   * @see K_means_center, K_means_clustering
   */
  class K_means_point
  {
  public:
    double data;      /**< Location of the point */
    int    center_id; /**< Closest center to the point */
    K_means_point(double data, int center_id = -1) : data(data),
      center_id(center_id) {
    }
    /**
     * Finds closest center and adds itself to the closest center's mean.
     * @param centers available centers
     * @return true if #center_id is changed (i.e. new center != previous center)
     */
    bool calculate_new_center(std::vector<K_means_center>& centers) {
      int new_center_id = 0;
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
  };

public:
  /** Types of algorithms for random center selection. */
  enum Initialization_types {
    RANDOM_INITIALIZATION, /**< place initial centers randomly */
    PLUS_INITIALIZATION    /**< place initial centers using k-means++ algorithm */
  };

private:
  std::vector<K_means_center> centers;
  std::vector<K_means_point>  points;
  int  maximum_iteration;

  Initialization_types init_type;

public:
  /**
   * Constructs structures and runs the algorithm.
   * K-means algorithm is repeated number_of_run times, and the result which has minimum within cluster error is kept.
   * @param number_of_centers
   * @param data
   * @param init_type initialization type for random center selection
   * @param number_of_run number of times to repeat k-means algorithm
   * @param maximum_iteration maximum allowed iteration in a single k-means algorithm call
   */
  K_means_clustering(int number_of_centers,
                     const std::vector<double>& data,
                     Initialization_types init_type = PLUS_INITIALIZATION,
                     int number_of_run = CGAL_DEFAULT_NUMBER_OF_RUN,
                     int maximum_iteration = CGAL_DEFAULT_MAXIMUM_ITERATION)
    :
    points(data.begin(), data.end()), maximum_iteration(maximum_iteration),
    init_type(init_type) {
    srand(CGAL_DEFAULT_SEED); //(static_cast<unsigned int>(time(NULL)))
    calculate_clustering_with_multiple_run(number_of_centers, number_of_run);
    sort(centers.begin(), centers.end());
  }

  /**
   * Fills data_center by the id of the closest center for each point.
   * @param[out] data_centers
   */
  void fill_with_center_ids(std::vector<int>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      data_centers.push_back(point_it->center_id);
    }
  }

private:
  /**
   * Initializes centers by choosing random points from data.
   * @param number_of_centers
   */
  void initiate_centers_randomly(int number_of_centers) {
    centers.clear();
    for(int i = 0; i < number_of_centers; ++i) {
      double initial_mean = points[rand() % points.size()].data;
      if(!make_center(initial_mean)) {
        --i;
      }
    }
  }

  /**
   * Initializes centers by using K-means++ algorithm.
   * Probability of a point to become a center is proportional to its squared distance to the closest center.
   * @param number_of_centers
   */
  void initiate_centers_plus_plus(int number_of_centers) {
    centers.clear();
    std::vector<double> distance_square_cumulative(points.size());
    std::vector<double> distance_square(points.size(),
                                        (std::numeric_limits<double>::max)());
    // distance_square stores squared distance to the closest center for each point.
    // say, distance_square -> [ 0.1, 0.2, 0.3, 0.4, ... ]
    // then corresponding distance_square_cumulative -> [ 0.1, 0.3, 0.6, 1, ...]
    double initial_mean = points[rand() % points.size()].data;
    make_center(initial_mean);

    for(int i = 1; i < number_of_centers; ++i) {
      double cumulative_distance_square = 0.0;
      for(std::size_t j = 0; j < points.size(); ++j) {
        double new_distance = std::pow(centers.back().mean - points[j].data, 2);
        if(new_distance < distance_square[j]) {
          distance_square[j] = new_distance;
        }
        cumulative_distance_square += distance_square[j];
        distance_square_cumulative[j] = cumulative_distance_square;
      }
      double zero_one = rand() / (RAND_MAX + 1.0); // [0,1) random number
      double random_ds =  zero_one * (distance_square_cumulative.back());
      int selection_index = std::upper_bound(distance_square_cumulative.begin(),
                                             distance_square_cumulative.end(), random_ds)
                            - distance_square_cumulative.begin();
      double initial_mean = points[selection_index].data;
      if(!make_center(initial_mean)) {
        --i;
      }
    }
  }

  /**
   * Checks whether the parameter center is previosly included in the center list.
   * @param center to be checked against existent centers
   * @return true if there is any center in the center list which has same mean with the parameter center
   */
  bool is_already_center(const K_means_center& center) const {
    for(std::vector<K_means_center>::const_iterator it = centers.begin();
        it != centers.end(); ++it) {
      if(it->mean == center.mean) {
        return true;
      }
    }
    return false;
  }

  /**
   * Adds @a value as a center if it is not previously added
   * @param value to become a center
   * @return true if center addition is succesful
   */
  bool make_center(double value) {
    K_means_center new_center(value);
    if(is_already_center(new_center)) {
      return false;
    }
    centers.push_back(new_center);
    return true;
  }
  /**
   * Main entry point for k-means algorithm.
   * Iterates until convergence occurs (i.e. no point changes its center) or maximum iteration limit is reached.
   */
  void calculate_clustering() {
    int iteration_count = 0;
    bool any_center_changed = true;
    while(any_center_changed && iteration_count++ < maximum_iteration) {
      any_center_changed = false;
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
    }
  }

  /**
   * Calls calculate_clustering() @a number_of_run times,
   * and keeps the result which has minimum within cluster error.
   * @param number_of_centers
   * @param number_of_run
   * @see calculate_clustering(), within_cluster_sum_of_squares()
   */
  void calculate_clustering_with_multiple_run(int number_of_centers,
      int number_of_run) {
    std::vector<K_means_center> min_centers;
    double error = (std::numeric_limits<double>::max)();
    while(number_of_run-- > 0) {
      clear_center_ids();
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
    // By saving points (min_points) also, we can get rid of this part
    // But since centers are already converged this step will require one iteration
    // If it stopped since maximum_iteration is reached then this step will take some time
    centers = min_centers;
    calculate_clustering();
  }

  /**
   * Sum of squared distances between each point and the closest center to it.
   */
  double within_cluster_sum_of_squares() const {
    double sum = 0;
    for(std::vector<K_means_point>::const_iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      sum += std::pow(centers[point_it->center_id].mean - point_it->data, 2);
    }
    return sum;
  }

  void clear_center_ids() {
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      point_it->center_id = -1;
    }
  }
};
}//namespace internal
}//namespace CGAL
#undef CGAL_DEFAULT_SEED
#undef CGAL_DEFAULT_MAXIMUM_ITERATION
#undef CGAL_DEFAULT_NUMBER_OF_RUN

#endif //CGAL_SURFACE_MESH_SEGMENTATION_K_MEANS_CLUSTERING_H
