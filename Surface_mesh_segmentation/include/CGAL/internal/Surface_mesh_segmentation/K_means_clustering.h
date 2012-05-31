#ifndef CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H
#define CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H

#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>

namespace CGAL
{

class K_means_point;

class K_means_center
{
public:
  double mean;
  int    number_of_points;
protected:
  double new_mean;

public:
  K_means_center(double mean): mean(mean), new_mean(0), number_of_points(0) {
  }
  void add_point(double data) {
    ++number_of_points;
    new_mean += data;
  }
  void calculate_mean() {
    mean = new_mean / number_of_points;
    new_mean = 0;
    number_of_points = 0;
  }
  bool operator < (const K_means_center& center) {
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
  /* Finds closest center to point,
  /* Adds itself to the closest center's mean,
  /* Returns true if center is changed.*/
  bool calculate_new_center(std::vector<K_means_center>& centers) {
    int new_center_id = 0;
    double min_distance = fabs(centers[0].mean - data);
    for(int i = 1; i < centers.size(); ++i) {
      double new_distance = fabs(centers[i].mean - data);
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
  K_means_clustering(int number_of_centers, const std::vector<double>& data,
                     int number_of_run = 1000, int maximum_iteration = 500)
    : points(data.begin(), data.end()), maximum_iteration(maximum_iteration),
      is_converged(false) {
    calculate_clustering_with_multiple_run(number_of_centers, number_of_run);
  }

  void fill_with_center_ids(std::vector<int>& data_centers) {
    for(std::vector<K_means_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      data_centers.push_back(point_it->center_id);
    }
  }

  void initiate_centers(int number_of_centers, bool randomly = false) {
    centers.clear();
    if(!randomly) {
      for(int i = 0; i < number_of_centers; ++i) {
        double initial_mean = (i + 1.0) / (number_of_centers + 1.0);
        centers.push_back(K_means_center(initial_mean));
      }
    } else {
      srand(time(NULL));
      for(int i = 0; i < number_of_centers; ++i) {
        double initial_mean = points[rand() % points.size()].data;
        centers.push_back(K_means_center(initial_mean));
      }
    }
    sort(centers.begin(), centers.end());
  }

  void calculate_clustering() {
    int iteration_count = 0;
    is_converged = false;
    do {
      for(std::vector<K_means_point>::iterator point_it = points.begin();
          point_it != points.end(); ++point_it) {
        bool center_changed = point_it->calculate_new_center(centers);
        is_converged |= center_changed;
      }
      for(std::vector<K_means_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it) {
        center_it->calculate_mean();
      }
    } while(!is_converged && ++iteration_count < maximum_iteration);
  }

  void calculate_clustering_with_multiple_run(int number_of_centers,
      int number_of_run) {
    initiate_centers(number_of_centers);
    calculate_clustering();

    std::vector<K_means_center> min_centers = centers;
    double error = sum_of_squares();
    while(--number_of_run > 0) {
      initiate_centers(number_of_centers, true);
      /* here, clearing center_ids of points might be neccessary */
      calculate_clustering();
      double new_error = sum_of_squares();
      if(error > new_error) {
        error = new_error;
        min_centers = centers;
      }
    }
    /* By saving points (min_points) also, we can get rid of this part */
    /* But since centers are already converged this step will require one iteration */
    centers = min_centers;
    calculate_clustering();
  }
  double sum_of_squares() const {
    double sum = 0;
    for(std::vector<K_means_center>::const_iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      for(std::vector<K_means_point>::const_iterator point_it = points.begin();
          point_it != points.end(); ++point_it) {
        sum += pow(center_it->mean - point_it->data, 2);
      }
    }
    return sum;
  }
};
}//namespace CGAL
#endif //CGAL_SEGMENTATION_K_MEANS_CLUSTERING_H
