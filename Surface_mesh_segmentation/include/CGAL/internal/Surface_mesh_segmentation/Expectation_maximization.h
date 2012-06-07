#ifndef CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
#define CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
/* NEED TO BE DONE */
/* About implementation:
 * Calculating probability multiple times, solution: storing matrix(cluster, point) for probability values.
 * About safe division: where centers have no members, in graph cut it may happen
 */
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <limits>

#include <fstream>
#include <iostream>
//#define ACTIVATE_SEGMENTATION_EM_DEBUG
#ifdef ACTIVATE_SEGMENTATION_EM_DEBUG
#define SEG_DEBUG(x) x;
#else
#define SEG_DEBUG(x)
#endif

namespace CGAL
{

class Gaussian_point;

class Gaussian_center
{
public:
  double mean;
  double deviation;
  double mixing_coefficient;

  Gaussian_center(): mean(0), deviation(0), mixing_coefficient(1.0) {
  }
  Gaussian_center(double mean, double deviation, double mixing_coefficient)
    : mean(mean), deviation(deviation), mixing_coefficient(mixing_coefficient) {
  }
  double probability(double x) const {
    double e_over = -0.5 * pow((x - mean) / deviation, 2);
    return exp(e_over) / deviation;
  }
  double probability_proportional(double x) const {
    double e_over = -0.5 * pow((x - mean) / deviation, 2);
    return exp(e_over);
  }
  bool operator < (const Gaussian_center& center) const {
    return mean < center.mean;
  }
  void calculate_parameters(const std::vector<Gaussian_point>& points);
};

class Gaussian_point
{
public:
  double data;
  double total_membership;

  Gaussian_point(double data): data(data), total_membership(0.0) {
  }
  void calculate_total_membership(const std::vector<Gaussian_center>& centers) {
    total_membership = 0.0;
    for(std::vector<Gaussian_center>::const_iterator it = centers.begin();
        it != centers.end(); ++it) {
      total_membership += it->probability(data) * it->mixing_coefficient;
    }
  }
};

inline void Gaussian_center::calculate_parameters(const
    std::vector<Gaussian_point>& points)
{
  /* Calculate new mean */
  double new_mean = 0.0, total_membership = 0.0;
  for(std::vector<Gaussian_point>::const_iterator it = points.begin();
      it != points.end(); ++it) {
    double membership = (probability(it->data) * mixing_coefficient) /
                        it->total_membership;
    new_mean += membership * it->data;
    total_membership += membership;
  }
  new_mean /= total_membership;
  /* Calculate new deviation */
  double new_deviation = 0.0;
  for(std::vector<Gaussian_point>::const_iterator it = points.begin();
      it != points.end(); ++it) {
    double membership = (probability(it->data) * mixing_coefficient) /
                        it->total_membership;
    new_deviation += membership * pow(it->data - new_mean, 2);
  }
  new_deviation = sqrt(new_deviation/total_membership);
  /* Calculate new mixing coefficient */
  mixing_coefficient = total_membership;
  deviation = new_deviation;
  mean = new_mean;
}

class Expectation_maximization
{
public:
  std::vector<Gaussian_center> centers;
  std::vector<Gaussian_point>  points;
  double threshold;
  int  maximum_iteration;
  bool is_converged;
protected:
  int seed;

public:
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           const std::vector<int>& initial_centers = std::vector<int>(),
                           int number_of_runs = 50, int maximum_iteration = 100)
    : points(data.begin(), data.end()), threshold(1e-4),
      maximum_iteration(maximum_iteration), is_converged(false), seed(time(NULL)) {
    srand(seed);
    if(initial_centers.empty()) {
      calculate_fitting_with_multiple_run(number_of_centers, number_of_runs);
    } else {
      initiate_centers(number_of_centers, initial_centers);
      calculate_fitting();
    }
  }

  void fill_with_center_ids(std::vector<int>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<Gaussian_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      double max_likelihood = 0.0;
      int max_center = 0, center_counter = 0;
      for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it, center_counter++) {
        double likelihood = center_it->mixing_coefficient * center_it->probability(
                              point_it->data);
        if(max_likelihood < likelihood) {
          max_likelihood = likelihood;
          max_center = center_counter;
        }
      }
      data_centers.push_back(max_center);
    }
  }
protected:
  void initiate_centers_randomly(int number_of_centers) {
    centers.clear();
    /* Randomly generate means of centers */
    //vector<int> center_indexes;
    double initial_mixing_coefficient = 1.0;
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    for(int i = 0; i < number_of_centers; ++i) {
      double random_index = rand() % points.size();
      //center_indexes.push_back(random_index);
      Gaussian_point mean_point = points[random_index];
      double initial_mean = mean_point.data;
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
    //write_random_centers("center_indexes.txt", center_indexes);
  }

  void initiate_centers_uniformly(int number_of_centers) {
    centers.clear();
    /* Uniformly generate centers */
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    double initial_mixing_coefficient = 1.0;
    for(int i = 0; i < number_of_centers; ++i) {
      double initial_mean = (i + 1.0) / (number_of_centers + 1.0);
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
  }

  void initiate_centers_from_memberships(int number_of_centers,
                                         const std::vector<int>& initial_centers) {
    centers.clear();
    /* Calculate mean */
    int number_of_point = initial_centers.size();
    centers = std::vector<Gaussian_center>(number_of_centers);
    std::vector<int> member_count(number_of_centers, 0);

    for(int i = 0; i < number_of_point; ++i) {
      int center_id = initial_centers[i];
      double data = points[i].data;
      centers[center_id].mean += data;
      member_count[center_id] += 1;
    }
    /* Assign mean, and mixing coef */
    for(int i = 0; i < number_of_centers; ++i) {
      centers[i].mean /= member_count[i];
      centers[i].mixing_coefficient =  member_count[i] / static_cast<double>
                                       (number_of_point);
    }
    /* Calculate deviation */
    for(int i = 0; i < number_of_point; ++i) {
      int center_id = initial_centers[i];
      double data = points[i].data;
      centers[center_id].deviation += pow(data - centers[center_id].mean, 2);
    }
    for(int i = 0; i < number_of_centers; ++i) {
      centers[i].deviation = sqrt(centers[i].deviation / member_count[i]);
    }
    sort(centers.begin(), centers.end());
  }

  void initiate_centers(int number_of_centers,
                        const std::vector<int>& initial_centers) {
    if(initial_centers.empty()) {
      initiate_centers_uniformly(number_of_centers);
    } else {
      initiate_centers_from_memberships(number_of_centers, initial_centers);
    }
  }
  /*Calculates total membership values for a point */
  void calculate_membership() {
    for(std::vector<Gaussian_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      point_it->calculate_total_membership(centers);
    }
  }
  double calculate_deviation(const Gaussian_point& point) const {
    double deviation = 0.0;
    for(std::vector<Gaussian_point>::const_iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      deviation += pow(point_it->data - point.data, 2);
    }
    return sqrt(deviation / points.size());
  }
  /*Calculates new parameter values for each cluster */
  void calculate_parameters() {
    for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      center_it->calculate_parameters(points);
    }
  }
  /*Calculates how much this adjustment is likely to represent data*/
  double calculate_likelihood() {
    double likelihood = 0.0;
    for(std::vector<Gaussian_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      double point_likelihood = 0.0;
      for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it) {
        point_likelihood += center_it->mixing_coefficient * center_it->probability(
                              point_it->data);
      }
      likelihood += log(point_likelihood);
    }
    return likelihood;
  }

  double iterate() {
    calculate_membership();
    calculate_parameters();
    return calculate_likelihood();
  }

  double calculate_fitting() {
    double likelihood = -(std::numeric_limits<double>::max)(), prev_likelihood;
    int iteration_count = 0;
    is_converged = false;
    while(!is_converged && iteration_count++ < maximum_iteration) {
      prev_likelihood = likelihood;
      likelihood = iterate();
      is_converged = likelihood - prev_likelihood < threshold * fabs(likelihood);
    }
    SEG_DEBUG(std::cout << likelihood << " " << iteration_count << std::endl)
    return likelihood;
  }

  void calculate_fitting_with_multiple_run(int number_of_centers,
      int number_of_run) {
    double max_likelihood = -(std::numeric_limits<double>::max)();
    std::vector<Gaussian_center> max_centers;

    while(number_of_run-- > 0) {
      initiate_centers_randomly(number_of_centers);

      double likelihood = calculate_fitting();
      write_center_parameters("center_paramters.txt");
      if(likelihood > max_likelihood) {
        max_centers = centers;
        max_likelihood = likelihood;
      }
    }
    centers = max_centers;
  }

  void write_random_centers(const char* file_name,
                            const std::vector<int>& center_indexes) {
    std::ofstream output(file_name, std::ios_base::app);
    for(std::vector<int>::const_iterator it = center_indexes.begin();
        it != center_indexes.end(); ++it) {
      output << points[(*it)].data << std::endl;
    }
    output.close();
  }

  void write_center_parameters(const char* file_name) {
    std::ofstream output(file_name, std::ios_base::app);
    for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      output << "mean: " << center_it->mean << " dev: " << center_it->deviation
             << " mix: " << center_it->mixing_coefficient << std::endl;
    }
  }
};
}//namespace CGAL
#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H