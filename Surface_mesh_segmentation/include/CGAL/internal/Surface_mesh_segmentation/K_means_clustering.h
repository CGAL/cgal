#ifndef CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H
#define CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H

#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <algorithm>

#define CGAL_DEFAULT_MAXIMUM_ITERATION 15
#define CGAL_DEFAULT_NUMBER_OF_RUN 20
#define CGAL_DEFAULT_SEED 1340818006

//#define ACTIVATE_SEGMENTATION_K_MEANS_DEBUG
#ifdef ACTIVATE_SEGMENTATION_K_MEANS_DEBUG
#define SEG_DEBUG(x) x;
#else
#define SEG_DEBUG(x)
#endif

namespace CGAL
{
namespace internal
{

class K_means_center
{
public:
  double mean;
protected:
  double new_mean;
  int    new_number_of_points;

public:
  K_means_center(double mean): mean(mean), new_mean(0), new_number_of_points(0) {
  }
  void add_point(double data) {
    ++new_number_of_points;
    new_mean += data;
  }
  void calculate_mean() {
    mean = new_mean / new_number_of_points;
    new_number_of_points = new_mean = 0;
  }
  bool operator < (const K_means_center& center) const {
    return mean < center.mean;
  }
};

class K_means_point
{
public:
  double data;
  int    center_id;
  K_means_point(double data, int center_id = -1) : data(data),
    center_id(center_id) {
  }
  /**
   * Finds closest center,
   * Adds itself to the closest center's mean,
   * Returns true if center is changed.
   */
  bool calculate_new_center(std::vector<K_means_center>& centers) {
    int new_center_id = 0;
    double min_distance = std::fabs(centers[0].mean - data);
    for(std::size_t i = 1; i < centers.size(); ++i) {
      double new_distance = std::fabs(centers[i].mean - data);
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

class K_means_clustering
{
public:
  std::vector<K_means_center> centers;
  std::vector<K_means_point>  points;
  int  maximum_iteration;
  bool is_converged;
protected:
  unsigned int seed;

public:
  K_means_clustering(int number_of_centers,
                     const std::vector<double>& data,
                     int number_of_run = CGAL_DEFAULT_NUMBER_OF_RUN,
                     int maximum_iteration = CGAL_DEFAULT_MAXIMUM_ITERATION)
    :
    points(data.begin(), data.end()), maximum_iteration(maximum_iteration),
    is_converged(false),
    seed(static_cast<unsigned int>(time(NULL))) {
#if 0
    srand(seed);
#else
    srand(CGAL_DEFAULT_SEED);
#endif
    calculate_clustering_with_multiple_run(number_of_centers, number_of_run);
  }
  /**
   * For each data point, data_center is filled by its center's id.
   */
  void fill_with_center_ids(std::vector<int>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      data_centers.push_back(point_it->center_id);
    }
  }

  void initiate_centers_uniformly(int number_of_centers) {
    centers.clear();
    for(int i = 0; i < number_of_centers; ++i) {
      double initial_mean = (i + 1.0) / (number_of_centers + 1.0);
      centers.push_back(K_means_center(initial_mean));
    }
    sort(centers.begin(), centers.end());
  }
  /* Forgy method */
  void initiate_centers_randomly(int number_of_centers) {
    centers.clear();
    for(int i = 0; i < number_of_centers; ++i) {
      double initial_mean = points[rand() % points.size()].data;
      K_means_center new_center(initial_mean);
      if(is_already_center(new_center)) {
        --i;
      } else                              {
        centers.push_back(new_center);
      }
    }
    sort(centers.begin(), centers.end());
  }

  void initiate_centers_plus_plus(int number_of_centers) {
    centers.clear();
    std::vector<double> distance_square_cumulative(points.size());
    std::vector<double> distance_square(points.size(),
                                        (std::numeric_limits<double>::max)());
    // distance_square stores squared distance to closest center for each point.
    // say, distance_square -> [ 0.1, 0.2, 0.3, 0.4, ... ]
    // then corresponding distance_square_cumulative -> [ 0.1, 0.3, 0.6, 1, ...]
    double initial_mean = points[rand() % points.size()].data;
    centers.push_back(K_means_center(initial_mean));

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

      double random_ds = (rand() / static_cast<double>(RAND_MAX)) *
                         (distance_square_cumulative.back());
      int selection_index = lower_bound(distance_square_cumulative.begin(),
                                        distance_square_cumulative.end(), random_ds)
                            - distance_square_cumulative.begin();
      double initial_mean = points[selection_index].data;
      K_means_center new_center(initial_mean);
      if(is_already_center(new_center)) {
        --i;
      } else                              {
        centers.push_back(new_center);
      }
    }
    sort(centers.begin(), centers.end());
  }

  bool is_already_center(const K_means_center& center) const {
    for(std::vector<K_means_center>::const_iterator it = centers.begin();
        it != centers.end(); ++it) {
      if(it->mean == center.mean) {
        return true;
      }
    }
    return false;
  }

  void calculate_clustering() {
    int iteration_count = 0;
    bool any_center_changed = true;
    while(any_center_changed && iteration_count++ < maximum_iteration) {
      any_center_changed = false;
      /* For each point, calculate its new center */
      for(std::vector<K_means_point>::iterator point_it = points.begin();
          point_it != points.end(); ++point_it) {
        bool center_changed = point_it->calculate_new_center(centers);
        any_center_changed |= center_changed;
      }
      /* For each center, calculate its new mean */
      for(std::vector<K_means_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it) {
        center_it->calculate_mean();
      }
    }
    is_converged = !any_center_changed;
    SEG_DEBUG(std::cout << iteration_count << " " << (is_converged ? "converged" :
              "not converged") << std::endl)
  }

  void calculate_clustering_with_multiple_run(int number_of_centers,
      int number_of_run) {
    std::vector<K_means_center> min_centers;
    double error = (std::numeric_limits<double>::max)();
    while(number_of_run-- > 0) {
      clear_center_ids();
#if 0
      initiate_centers_randomly(number_of_centers);
#else
      initiate_centers_plus_plus(number_of_centers);
#endif

      calculate_clustering();
#if 0
      double new_error = sum_of_squares();
#else
      double new_error = within_cluster_sum_of_squares();
#endif
      if(error > new_error) {
        error = new_error;
        min_centers = centers;
      }
    }
    /* By saving points (min_points) also, we can get rid of this part */
    /* But since centers are already converged this step will require one iteration */
    /* If it stopped since maximum_iteration is reached then this step will take some time */
    centers = min_centers;
    calculate_clustering();
  }

  double sum_of_squares() const {
    double sum = 0;
    for(std::vector<K_means_center>::const_iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      for(std::vector<K_means_point>::const_iterator point_it = points.begin();
          point_it != points.end(); ++point_it) {
        sum += std::pow(center_it->mean - point_it->data, 2);
      }
    }
    return sum;
  }

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

#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H
