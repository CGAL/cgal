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

#define DEF_MAX_ITER  10
#define DEF_THRESHOLD 1e-3
#define USE_MATRIX    true
#define DEF_SEED 1340818006

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
  double log_probability_with_coef(double x) const {
    return  log(mixing_coefficient) - log(deviation)
            -0.5 * pow((x - mean) / deviation, 2);
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
  double final_likelihood;
protected:
  unsigned int seed;
  std::vector<std::vector<double> > membership_matrix;

public:
  // these two constructors will be removed!
  Expectation_maximization() { }
  Expectation_maximization(const std::vector<double>& data): points(data) { }
  /* For uniform initialization, with one run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           double threshold = DEF_THRESHOLD,  int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false), seed(0),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))),
      final_likelihood(-(std::numeric_limits<double>::max)()) {
    initiate_centers_uniformly(number_of_centers);
    fit();
  }
  /* For initialization with provided center ids per point, with one run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           std::vector<int>& initial_center_ids,
                           double threshold = DEF_THRESHOLD, int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false), seed(0),
      final_likelihood(-(std::numeric_limits<double>::max)()) {
    number_of_centers = initiate_centers_from_memberships(number_of_centers,
                        initial_center_ids);
    membership_matrix = std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()));
    fit();
    //write_center_parameters("centers_param.txt");
  }
  /* For initialization with random center selection (Forgy), with multiple run */
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           int number_of_runs,
                           double threshold = DEF_THRESHOLD,  int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false),
      seed(static_cast<unsigned int>(time(NULL))),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))),
      final_likelihood(-(std::numeric_limits<double>::max)()) {
#if 0
    srand(seed);
#else
    srand(DEF_SEED);
#endif
    fit_with_multiple_run(number_of_centers, number_of_runs);
    //write_center_parameters("centers_param.txt");
  }

  // Experimental!
  Expectation_maximization(int number_of_centers, const std::vector<double>& data,
                           int number_of_runs,
                           std::vector<std::vector<double> >& probabilities,
                           std::vector<int>& data_centers,
                           double threshold = DEF_THRESHOLD,  int maximum_iteration = DEF_MAX_ITER)
    : points(data), threshold(threshold), maximum_iteration(maximum_iteration),
      is_converged(false),
      seed(static_cast<unsigned int>(time(NULL))),
      membership_matrix(std::vector<std::vector<double> >(number_of_centers,
                        std::vector<double>(points.size()))),
      final_likelihood(-(std::numeric_limits<double>::max)()) {
#if 0
    srand(seed);
#else
    srand(DEF_SEED);
#endif
    average_with_multiple_run(number_of_centers, number_of_runs, probabilities,
                              data_centers);
  }

  void fill_with_center_ids(std::vector<int>& data_centers) {
    data_centers.reserve(points.size());
    for(std::vector<double>::iterator point_it = points.begin();
        point_it != points.end(); ++point_it) {
      double max_likelihood = 0.0;
      int max_center = 0, center_counter = 0;
      for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
          center_it != centers.end(); ++center_it, center_counter++) {
        double likelihood = center_it->probability_with_coef(*point_it);
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
    double epsilon = 1e-8;
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
        double probability = probabilities[center_i][point_i] / total_probability;
        probability += epsilon;
        probability = (std::min)(probability, 1.0);
        probability = -log(probability);
        probabilities[center_i][point_i] = (std::max)(probability,
                                           std::numeric_limits<double>::epsilon());
      }
    }
  }

  // experimental, going to be removed most prob !
  void refresh_parameters_and_probabilities(int number_of_centers,
      std::vector<int>& initial_center_ids,
      std::vector<std::vector<double> >& probabilities) {
    std::vector<bool> cluster_exist(number_of_centers, false);
    std::vector<int>  cluster_shift(number_of_centers, 0);
    for(std::vector<int>::iterator id_it = initial_center_ids.begin();
        id_it != initial_center_ids.end(); ++id_it) {
      cluster_exist[*id_it] = true;
    }

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
      if(!cluster_exist[i]) {
        continue;
      }
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
      if(!cluster_exist[i]) {
        continue;
      }
      centers[i].deviation = sqrt(centers[i].deviation / member_count[i]);
    }

    double epsilon = 1e
                     -5; // this epsilon should be consistent with epsilon in calculate_dihedral_angle_of_edge!
    probabilities = std::vector<std::vector<double> >
                    (number_of_centers, std::vector<double>(points.size()));
    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      double total_probability = 0.0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        if(!cluster_exist[center_i]) {
          probabilities[center_i][point_i] = 0;
          continue;
        }
        double probability = centers[center_i].probability(points[point_i]);
        total_probability += probability;
        probabilities[center_i][point_i] = probability;
      }
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double probability = probabilities[center_i][point_i] / total_probability;
        probability = (std::max)(probability, epsilon);
        probability = -log(probability);
        probabilities[center_i][point_i] = (std::max)(probability,
                                           1e-5); // this is another epsilon, edge can not hold 0 to source in graph-cut.
      }
    }
    //sort(centers.begin(), centers.end());
  }
protected:

  // Initialization functions for centers.

  void initiate_centers_randomly(int number_of_centers) {
    centers.clear();
    /* Randomly generate means of centers */
    double initial_mixing_coefficient = 1.0 / number_of_centers;
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    for(int i = 0; i < number_of_centers; ++i) {
      int random_index = rand() % points.size();
      double initial_mean = points[random_index];
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
    //write_random_centers("center_indexes.txt");
  }

  void initiate_centers_uniformly(int number_of_centers) {
    centers.clear();
    /* Uniformly generate means of centers */
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    double initial_mixing_coefficient = 1.0 / number_of_centers;
    for(int i = 0; i < number_of_centers; ++i) {
      double initial_mean = (i + 1.0) / (number_of_centers + 1.0);
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
  }

  void initiate_centers_plus_plus(int number_of_centers) {
    centers.clear();
    double initial_deviation = 1.0  / (2.0 * number_of_centers);
    double initial_mixing_coefficient = 1.0 / number_of_centers;

    std::vector<double> distance_square_cumulative(points.size());
    std::vector<double> distance_square(points.size(),
                                        (std::numeric_limits<double>::max)());
    double initial_mean = points[rand() % points.size()];
    centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                      initial_mixing_coefficient));

    for(int i = 1; i < number_of_centers; ++i) {
      double cumulative_distance_square = 0.0;
      for(std::size_t j = 0; j < points.size(); ++j) {
        double new_distance = pow(centers.back().mean - points[j], 2);
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
      double initial_mean = points[selection_index];
      centers.push_back(Gaussian_center(initial_mean, initial_deviation,
                                        initial_mixing_coefficient));
    }
    sort(centers.begin(), centers.end());
  }

  int initiate_centers_from_memberships(int number_of_centers,
                                        std::vector<int>& initial_center_ids) {
    /* For handling clusters that have 0 members (in initial_center_ids) */
    /* remove those clusters, and shift cluster ids. */
    std::vector<bool> cluster_exist(number_of_centers, false);
    std::vector<int>  cluster_shift(number_of_centers, 0);
    for(std::vector<int>::iterator id_it = initial_center_ids.begin();
        id_it != initial_center_ids.end(); ++id_it) {
      cluster_exist[*id_it] = true;
    }
    int shift = 0;
    for(int i = 0; i < number_of_centers; ++i) {
      if(!cluster_exist[i]) {
        --shift;
      }
      cluster_shift[i] = shift;
    }
    number_of_centers += shift;
    for(std::vector<int>::iterator id_it = initial_center_ids.begin();
        id_it != initial_center_ids.end(); ++id_it) {
      *id_it += cluster_shift[*id_it];
    }

    /* Calculate mean */
    int number_of_points = initial_center_ids.size();
    centers = std::vector<Gaussian_center>(number_of_centers);
    std::vector<int> member_count(number_of_centers, 0);

    for(int i = 0; i < number_of_points; ++i) {
      int center_id = initial_center_ids[i];
      centers[center_id].mean += points[i];
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
      centers[center_id].deviation += pow(points[i] - centers[center_id].mean, 2);
    }
    for(int i = 0; i < number_of_centers; ++i) {
      centers[i].deviation = sqrt(centers[i].deviation / member_count[i]);
    }
    sort(centers.begin(), centers.end());
    return number_of_centers;
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
      centers[center_i].mixing_coefficient = total_membership / points.size();
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

  // Both for E step, and likelihood step
  double calculate_likelihood_with_log_sum_exp() {
    /**
     * Calculate Log-likelihood
     * The trick (merely a trick) is while calculating log-likelihood, we also store probability results into matrix,
     * so that in next iteration we do not have to calculate them again.
     */
    double likelihood = 0.0;
    std::vector<double> max_probability(points.size(),
                                        -(std::numeric_limits<double>::max)());
    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double membership = centers[center_i].log_probability_with_coef(
                              points[point_i]);
        if(max_probability[point_i] < membership) {
          max_probability[point_i] = membership;
        }
        membership_matrix[center_i][point_i] = membership;
      }
      double total_membership = 0.0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        // here no need to multiply with exp(max_probability[point_i]) since it is going to be vanish
        // when we divide every element to 'total_membership'.
        double membership = exp(membership_matrix[center_i][point_i] -
                                max_probability[point_i]);
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
      double progress = likelihood - prev_likelihood;
      is_converged = progress > 0.0 && progress < threshold * std::fabs(likelihood);
    }
    SEG_DEBUG(std::cout << "likelihood: " <<  likelihood << "iteration: " <<
              iteration_count << std::endl)
    if(final_likelihood < likelihood) {
      final_likelihood = likelihood;
    }
    return likelihood;
  }

  // Calls fit() number_of_run times, and stores best found centers
  void fit_with_multiple_run(int number_of_centers, int number_of_run) {
    std::vector<Gaussian_center> max_centers;

    while(number_of_run-- > 0) {
#if 0
      initiate_centers_randomly(number_of_centers);
#else
      initiate_centers_plus_plus(number_of_centers);
#endif
      double likelihood = fit();
      if(likelihood == final_likelihood) {
        max_centers = centers;
      }
    }
    SEG_DEBUG(std::cout << "max likelihood: " << final_likelihood << std::endl)
    centers = max_centers;
    //write_center_parameters("centers_param.txt");
  }

  // Experimental!
  void average_with_multiple_run(int number_of_centers, int number_of_run,
                                 std::vector<std::vector<double> >& probabilities,
                                 std::vector<int>& data_centers) {
    double total_likelihood = 0.0;
    int temp_run = number_of_run;
    while(number_of_run-- > 0) {
      initiate_centers_randomly(number_of_centers);
      double likelihood = fit();
      total_likelihood += likelihood;

      for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
        std::vector<double> temp(centers.size());
        double total_probability = 0.0;
        for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
          double probability = centers[center_i].probability(points[point_i]);
          total_probability += probability;
          temp[center_i] = probability * likelihood;
        }
        for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
          probabilities[center_i][point_i] += temp[center_i]  / total_probability;
        }
      }
    }

    for(std::size_t point_i = 0; point_i < points.size(); ++point_i) {
      double max_probability = 0.0;
      int max_center = 0;
      for(std::size_t center_i = 0; center_i < centers.size(); ++center_i) {
        double probability = probabilities[center_i][point_i] / total_likelihood;
        if(max_probability < probability) {
          max_probability = probability;
          max_center = center_i;
        }
        probability = (std::max)(probability, 1e-8);
        probability = -log(probability);
        probabilities[center_i][point_i] = (std::max)(probability, 1e-8);
      }
      data_centers.push_back(max_center);
    }
  }

  void write_random_centers(const char* file_name) {
    //std::ofstream output(file_name, std::ios_base::app);
    std::ofstream output(file_name);
    for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      output << center_it->mean << std::endl;
    }
    output.close();
  }

  void write_center_parameters(const char* file_name) {
    //std::ofstream output(file_name, std::ios_base::app);
    std::ofstream output(file_name);
    for(std::vector<Gaussian_center>::iterator center_it = centers.begin();
        center_it != centers.end(); ++center_it) {
      output << "mean: " << center_it->mean << " dev: " << center_it->deviation
             << " mix: " << center_it->mixing_coefficient << std::endl;
    }
  }
};
}//namespace internal
}//namespace CGAL
#undef DEF_SEED
#undef DEF_MAX_ITER
#undef DEF_THRESHOLD
#undef USE_MATRIX
#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SEGMENTATION_EXPECTATION_MAXIMIZATION_H