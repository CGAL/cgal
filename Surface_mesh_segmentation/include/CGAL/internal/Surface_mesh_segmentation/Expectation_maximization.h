#ifndef CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
#define CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
/* NEED TO BE DONE */
/* About implementation:
/* Calculating probability multiple times */
/* About safe division: where centers have no members, in graph cut it may happen  */
#include <vector>
#include <cmath>
#include <algorithm>
//#define DIV_INV_SQRT_2_PI 0.3989422804
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
  bool operator < (const Gaussian_center& center) {
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
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           const std::vector<int>& initial_centers = std::vector<int>(),
                           int maximum_iteration = 100)
    : points(data.begin(), data.end()), threshold(1e-4),
      maximum_iteration(maximum_iteration), is_converged(false) {
    initiate_centers(number_of_centers, initial_centers);
    calculate_fitting();
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
  void initiate_centers(int number_of_centers,
                        const std::vector<int>& initial_centers) {
    if(initial_centers.empty()) {
      /* Uniformly generate centers */
      double initial_deviation = 1.0  / (2.0 * number_of_centers);
      double initial_mixing_coefficient = 1.0 / number_of_centers;
      for(int i = 0; i < number_of_centers; ++i) {
        double initial_mean = (i + 1.0) / (number_of_centers + 1.0);
        centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                          initial_mixing_coefficient));
      }
    } else {
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
    }
    sort(centers.begin(), centers.end());
  }
  /*Calculates total membership values for a point */
  void calculate_membership() {
    for(std::vector<Gaussian_point>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      point_it->calculate_total_membership(centers);
    }
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
  void calculate_fitting() {
    double prev_likelihood = 0.0;
    double likelihood = iterate();
    int iteration_count = 0;
    is_converged = false;
    do {
      prev_likelihood = likelihood;
      likelihood = iterate();
      is_converged = likelihood - prev_likelihood < threshold * likelihood;
    } while(!is_converged && ++iteration_count < maximum_iteration);
  }
};
}//namespace CGAL
#endif //CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H