#ifndef CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
#define CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H
/* NEED TO BE DONE */
/* About implementation:
 * + There are a lot of parameters (with default values) and initialization types,
 * so I am planning to use whether 'Named Parameter Idiom' or Boost Parameter Library...
 *
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>

#define DEF_MAX_ITER  100
#define DEF_THRESHOLD 1e-4
#define USE_MATRIX    true

#define ACTIVATE_SEGMENTATION_EM_DEBUG
#ifdef ACTIVATE_SEGMENTATION_EM_DEBUG
#define SEG_DEBUG(x) x;
#else
#define SEG_DEBUG(x)
#endif

namespace CGAL
{
namespace internal
{

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
  double probability_with_coef(double x) const {
    return probability(x) * mixing_coefficient;
  }
  bool operator < (const Gaussian_center& center) const {
    return mean < center.mean;
  }
};


class Expectation_maximization
{
public:
  std::vector<Gaussian_center> centers;
  std::vector<double>  points;
  double threshold;
  int  maximum_iteration;
  bool is_converged;
protected:
  unsigned int seed;
  std::vector<std::vector<double> > membership_matrix;

public:
  /* For uniform initialization, with one run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           double threshold = DEF_THRESHOLD,  int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false), seed(0),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))) {
    initiate_centers_uniformly(number_of_centers);
    fit();
  }
  /* For initialization with provided center ids per point, with one run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           const std::vector<int>& initial_center_ids,
                           double threshold = DEF_THRESHOLD, int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false), seed(0),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))) {
    initiate_centers_from_memberships(number_of_centers, initial_center_ids);
    fit();
  }
  /* For initialization with random center selection (Forgy), with multiple run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           int number_of_runs,
                           double threshold = DEF_THRESHOLD,  int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false),
      seed(static_cast<unsigned int>(time(NULL))),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))) {
    srand(seed);
    fit_with_multiple_run(number_of_centers, number_of_runs);
  }

  void fill_with_center_ids(std::vector<int>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<double>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      double max_likelihood = 0.0;
      int max_center = 0, center_counter = 0;
      for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it, center_counter++) {
        double likelihood = center_it->probability(*point_it);
        if(max_likelihood < likelihood) {
          max_likelihood = likelihood;
          max_center = center_counter;
        }
      }
      data_centers.push_back(max_center);
    }
  }

  void fill_with_minus_log_probabilities(std::vector<std::vector<double> >&
                                         probabilities) {
    double epsilon = 1e
                     -8; // this epsilon should be consistent with epsilon in calculate_dihedral_angle_of_edge!
    probabilities = std::vector<std::vector<double> >
                    (centers.size(), std::vector<double>(points.size()));
    for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
      double sum = 0.0;

      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double probability = centers[center_i].probability(points[point_i]);
        sum += probability;
        probabilities[center_i][point_i] = probability;
      }
#if 0
      // pdf values scaled so that their sum will equal to 1.
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double probability = probabilities[center_i][point_i] / sum;
        probability = (CGAL::max)(probability, epsilon);
        probabilities[center_i][point_i] = -log(probability);
      }
#else
      //pdf values scaled between [0-1]
      std::pair<std::vector<double>::iterator, std::vector<double>::iterator>
      min_max_pair =
        CGAL::min_max_element(probabilities[center_i].begin(),
                              probabilities[center_i].end());
      double max_value = *min_max_pair.second, min_value = *min_max_pair.first;
      double max_min_dif = max_value - min_value;
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double probability = probabilities[center_i][point_i];
        probability = (probability - min_value) / max_min_dif;
        probability = (CGAL::max)(probability, epsilon);
        probabilities[center_i][point_i] = -log(probability);
      }
#endif

    }
  }
protected:

  // Initialization functions for centers.

  void initiate_centers_randomly(int number_of_centers) {
    centers.clear();
    /* Randomly generate means of centers */
    //vector<int> center_indexes;
    double initial_mixing_coefficient = 1.0;
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    for(int i = 0; i < number_of_centers; ++i) {
      int random_index = rand() % points.size();
      //center_indexes.push_back(random_index);
      double initial_mean = points[random_index];
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
    //write_random_centers("center_indexes.txt", center_indexes);
  }

  void initiate_centers_uniformly(int number_of_centers) {
    centers.clear();
    /* Uniformly generate means of centers */
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
                                         const std::vector<int>& initial_center_ids) {
    /* Calculate mean */
    int number_of_points = initial_center_ids.size();
    centers = std::vector<Gaussian_center>(number_of_centers);
    std::vector<int> member_count(number_of_centers, 0);

    for(int i = 0; i < number_of_points; ++i) {
      int center_id = initial_center_ids[i];
      double data = points[i];
      centers[center_id].mean += data;
      member_count[center_id] += 1;
    }
    /* Assign mean, and mixing coef */
    for(int i = 0; i < number_of_centers; ++i) {
      centers[i].mean /= member_count[i];
      centers[i].mixing_coefficient =  member_count[i] / static_cast<double>
                                       (number_of_points);
    }
    /* Calculate deviation */
    for(int i = 0; i < number_of_points; ++i) {
      int center_id = initial_center_ids[i];
      double data = points[i];
      centers[center_id].deviation += pow(data - centers[center_id].mean, 2);
    }
    for(int i = 0; i < number_of_centers; ++i) {
      centers[i].deviation = sqrt(centers[i].deviation / member_count[i]);
    }
    sort(centers.begin(), centers.end());
  }

  //Main steps of EM

  // M step
  void calculate_parameters() {
    for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
      // Calculate new mean
      double new_mean = 0.0, total_membership = 0.0;
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double membership = membership_matrix[center_i][point_i];
        new_mean += membership * points[point_i];
        total_membership += membership;
      }
      new_mean /= total_membership;
      // Calculate new deviation
      double new_deviation = 0.0;
      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        double membership = membership_matrix[center_i][point_i];
        new_deviation += membership * pow(points[point_i] - new_mean, 2);
      }
      new_deviation = sqrt(new_deviation/total_membership);
      // Assign new parameters
      centers[center_i].mixing_coefficient = total_membership;
      centers[center_i].deviation = new_deviation;
      centers[center_i].mean = new_mean;
    }
  }

  // Both for E step, and likelihood step
  double calculate_likelihood() {
    /**
     * Calculate Log-likelihood
     * The trick (merely a trick) is while calculating log-likelihood, we also store probability results into matrix,
     * so that in next iteration we do not have to calculate them again.
     */
    double likelihood = 0.0;
    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      double total_membership = 0.0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double membership = centers[center_i].probability_with_coef(points[point_i]);
        total_membership += membership;
        membership_matrix[center_i][point_i] = membership;
      }
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        membership_matrix[center_i][point_i] /= total_membership;
      }
      likelihood += log(total_membership);
    }
    return likelihood;
  }

  // One iteration of EM
  double iterate(bool first_iteration) {
    // E-step
    if(first_iteration) {
      calculate_likelihood();
    }
    // M-step
    calculate_parameters();
    // Likelihood step
    return calculate_likelihood();
  }

  // Fitting function, iterates until convergence occurs or maximum iteration limit is reached
  double fit() {
    double likelihood = -(std::numeric_limits<double>::max)(), prev_likelihood;
    int iteration_count = 0;
    is_converged = false;
    while(!is_converged && iteration_count++ < maximum_iteration) {
      prev_likelihood = likelihood;
      likelihood = iterate(iteration_count == 1);
      is_converged = likelihood - prev_likelihood < threshold * fabs(likelihood);
    }
    //SEG_DEBUG(std::cout << likelihood << " " << iteration_count << std::endl)
    return likelihood;
  }

  // Calls fit() number_of_run times, and stores best found centers
  void fit_with_multiple_run(int number_of_centers, int number_of_run) {
    double max_likelihood = -(std::numeric_limits<double>::max)();
    std::vector<Gaussian_center> max_centers;

    while(number_of_run-- > 0) {
      initiate_centers_randomly(number_of_centers);

      double likelihood = fit();
      //write_center_parameters("center_paramters.txt");
      if(likelihood > max_likelihood) {
        max_centers = centers;
        max_likelihood = likelihood;
      }
    }
    SEG_DEBUG(std::cout << "max likelihood: " << max_likelihood << std::endl)
    centers = max_centers;
  }

  void write_random_centers(const char* file_name,
                            const std::vector<int>& center_indexes) {
    std::ofstream output(file_name, std::ios_base::app);
    for(std::vector<int>::const_iterator it = center_indexes.begin();
        it != center_indexes.end(); ++it) {
      output << points[(*it)] << std::endl;
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
}//namespace internal
}//namespace CGAL
#undef DEF_MAX_ITER
#undef DEF_THRESHOLD
#undef USE_MATRIX
#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H