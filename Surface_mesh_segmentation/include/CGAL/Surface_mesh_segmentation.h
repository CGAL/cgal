#ifndef CGAL_SURFACE_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_H
/* NEED TO BE DONE
 * About implementation:
 * +) I am not using BGL, as far as I checked there is a progress on BGL redesign
 *    (https://cgal.geometryfactory.com/CGAL/Members/wiki/Features/BGL) which introduces some features
 *    for face-based traversal / manipulation by FaceGraphs
 * +) Deciding on which parameters will be taken from user
 * +) Make it more readable: calculate_sdf_value_of_facet function.
 *
 * About paper (and correctness / efficiency etc.):
 * +) Weighting ray distances with inverse of their angles: not sure how to weight exactly
 * +) Anisotropic smoothing: have no idea what it is exactly, should read some material (google search is not enough)
 * +) Deciding how to generate rays in cone: for now using "polar angle" and "accept-reject (square)" and "concentric mapping" techniques
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

//#include "Expectation_maximization.h"
//#include "K_means_clustering.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h>
#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/utility.h>

#define LOG_5 1.60943791
#define NORMALIZATION_ALPHA 4.0
#define ANGLE_ST_DEV_DIVIDER 3.0

namespace CGAL
{

template <class Polyhedron>
class Surface_mesh_segmentation
{
//type definitions
public:
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Kernel::Vector_3   Vector;
  typedef typename Kernel::Point_3    Point;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle   Facet_handle;
protected:
  typedef typename Kernel::Ray_3 Ray;
  typedef typename CGAL::AABB_polyhedron_triangle_primitive<Kernel, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<CGAL::AABB_traits<Kernel, Primitive> >
  Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;

  typedef std::map<Facet_handle, double> Face_value_map;
  typedef std::map<Facet_handle, int>    Face_center_map;
  /*Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = angle with cone-normal (weight). */
  typedef CGAL::Triple<double, double, double> Disk_sample;
  typedef std::vector<CGAL::Triple<double, double, double> > Disk_samples_list;

  template <typename ValueTypeName>
  struct compare_pairs {
    bool operator()(const ValueTypeName& v1,const ValueTypeName& v2) const {
      return v1.second < v2.second;
    }
  };

  template <typename ValueTypeName>
  struct compare_pairs_using_first {
    bool operator()(const ValueTypeName& v1,const ValueTypeName& v2) const {
      return v1.first < v2.first;
    }
  };
//member variables
public:
  Polyhedron* mesh;

  Face_value_map  sdf_values;
  Face_center_map centers;

  double cone_angle;
  int    number_of_rays_sqrt;
  int    number_of_centers;
  /*Store sampled points from disk, just for once */
  Disk_samples_list disk_samples;

  std::ofstream log_file;

  //std::vector<int> center_memberships_temp;

//member functions
public:
  Surface_mesh_segmentation(Polyhedron* mesh,
                            int number_of_rays_sqrt = 9, double cone_angle = (2.0 / 3.0) * CGAL_PI,
                            int number_of_centers = 2);

  void calculate_sdf_values();

  double calculate_sdf_value_of_facet (const Facet_handle& facet,
                                       const Point& center,
                                       const Vector& normal_const, const Tree& tree) const;
  void cast_and_return_minimum (const Ray& ray, const Tree& tree,
                                const Facet_handle& facet,
                                bool& is_found, double& min_distance) const;

  double calculate_sdf_value_from_rays (std::vector<double>& ray_distances,
                                        std::vector<double>& ray_weights) const;
  double calculate_sdf_value_from_rays_with_mean (std::vector<double>&
      ray_distances, std::vector<double>& ray_weights) const;
  double calculate_sdf_value_from_rays_with_trimmed_mean (
    std::vector<double>& ray_distances, std::vector<double>& ray_weights) const;

  void disk_sampling_rejection();
  void disk_sampling_polar_mapping();
  void disk_sampling_concentric_mapping();

  void normalize_sdf_values();
  void smooth_sdf_values();

  void apply_GMM_fitting();
  void apply_K_means_clustering();
  void apply_GMM_fitting_with_K_means_init();

  void write_sdf_values(const char* file_name);
  void read_sdf_values(const char* file_name);

};

template <class Polyhedron>
inline Surface_mesh_segmentation<Polyhedron>::Surface_mesh_segmentation(
  Polyhedron* mesh, int number_of_rays_sqrt, double cone_angle,
  int number_of_centers)
  : mesh(mesh), cone_angle(cone_angle), number_of_rays_sqrt(number_of_rays_sqrt),
    number_of_centers(number_of_centers), log_file("log_file.txt")
{
  disk_sampling_concentric_mapping();
  calculate_sdf_values();
  apply_GMM_fitting_with_K_means_init();
  //write_sdf_values("sdf_values_sample_cactus.txt");
  //read_sdf_values("sdf_values_sample_cactus.txt");
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::calculate_sdf_values()
{
  sdf_values.clear();
  Tree tree(mesh->facets_begin(), mesh->facets_end());
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    CGAL_precondition(facet_it->is_triangle()); //Mesh should contain triangles.

    Point v1 = facet_it->halfedge()->vertex()->point();
    Point v2 = facet_it->halfedge()->next()->vertex()->point();
    Point v3 = facet_it->halfedge()->next()->next()->vertex()->point();
    Point center  = CGAL::centroid(v1, v2, v3);
    Vector normal = CGAL::unit_normal(v1, v2,
                                      v3) * -1.0; //Assuming triangles are CCW oriented.
    //SL: cone angle and number of rays should be parameters.
    double sdf = calculate_sdf_value_of_facet(facet_it, center, normal, tree);
    sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf));
  }
  normalize_sdf_values();
  smooth_sdf_values();
}

template <class Polyhedron>
inline double
Surface_mesh_segmentation<Polyhedron>::calculate_sdf_value_of_facet(
  const Facet_handle& facet, const Point& center,
  const Vector& normal_const, const Tree& tree) const
{
  typename Kernel::Plane_3 plane(center, normal_const);
  Vector v1 = plane.base1();
  Vector v2 = plane.base2();
  v1 = v1 / CGAL::sqrt(v1.squared_length());
  v2 = v2 / CGAL::sqrt(v2.squared_length());

  int ray_count = number_of_rays_sqrt * number_of_rays_sqrt;

  std::vector<double> ray_distances, ray_weights;
  ray_distances.reserve(ray_count);
  ray_weights.reserve(ray_count);

  double length_of_normal = 1.0 / tan(cone_angle / 2);
  Vector normal = normal_const * length_of_normal;

  for(Disk_samples_list::const_iterator sample_it = disk_samples.begin();
      sample_it != disk_samples.end(); ++sample_it) {
    Vector disk_vector = v1 * sample_it->first + v2 * sample_it->second;
    Ray ray(center, normal + disk_vector);
    double min_distance;
    bool is_found;
    cast_and_return_minimum(ray, tree, facet, is_found, min_distance);
    if(!is_found) {
      continue;
    }

    ray_weights.push_back(sample_it->third);
    ray_distances.push_back(min_distance);
  }
  return calculate_sdf_value_from_rays_with_trimmed_mean(ray_distances,
         ray_weights);
}

template <class Polyhedron>
void Surface_mesh_segmentation<Polyhedron>::cast_and_return_minimum(
  const Ray& ray, const Tree& tree, const Facet_handle& facet,
  bool& is_found, double& min_distance) const
{
  std::list<Object_and_primitive_id> intersections;
  tree.all_intersections(ray, std::back_inserter(intersections));
  Vector min_i_ray;
  typename Tree::Primitive_id min_id;
  is_found = false;
  for(typename std::list<Object_and_primitive_id>::iterator op_it =
        intersections.begin();
      op_it != intersections.end() ; ++op_it) {
    CGAL::Object object   = op_it->first;
    typename Tree::Primitive_id id = op_it->second;
    Point i_point;
    if(id == facet)                    {
      continue;  //Since center is located on related facet, we should skip it if there is an intersection with it.
    }
    if(!CGAL::assign(i_point, object)) {
      continue;  //What to do here (in case of intersection object is a segment), I am not sure ???
    }
    Vector i_ray = (ray.source() - i_point);
    double new_distance = CGAL::sqrt(i_ray.squared_length());
    if(!is_found || new_distance < min_distance) {
      min_distance = new_distance;
      min_id = id;
      min_i_ray = i_ray;
      is_found = true;
    }
  }
  if(!is_found) {
    return;
  }
  Point min_v1 = min_id->halfedge()->vertex()->point();
  Point min_v2 = min_id->halfedge()->next()->vertex()->point();
  Point min_v3 = min_id->halfedge()->next()->next()->vertex()->point();
  Vector min_normal = CGAL::normal(min_v1, min_v2, min_v3) * -1.0;

  if(CGAL::angle(CGAL::ORIGIN + min_i_ray, Point(CGAL::ORIGIN),
                 CGAL::ORIGIN + min_normal) != CGAL::ACUTE) {
    is_found = false;
  }
}

template <class Polyhedron>
inline double
Surface_mesh_segmentation<Polyhedron>::calculate_sdf_value_from_rays(
  std::vector<double>& ray_distances,
  std::vector<double>& ray_weights) const
{
  int accepted_ray_count = ray_distances.size();
  if(accepted_ray_count == 0)      {
    return 0.0;
  } else if(accepted_ray_count == 1) {
    return ray_distances[0];
  }
  /* local variables */
  double total_weights = 0.0, total_distance = 0.0;
  double median_sdf = 0.0, mean_sdf = 0.0, st_dev = 0.0;
  /* Calculate mean sdf */
  std::vector<double>::iterator w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end();
      ++dist_it, ++w_it) {
    total_distance += (*dist_it) * (*w_it);
    total_weights += (*w_it);
  }
  mean_sdf = total_distance / total_weights;
  total_weights = 0.0;
  total_distance = 0.0;
  /* Calculate median sdf */
  int half_ray_count = accepted_ray_count / 2;
  std::nth_element(ray_distances.begin(), ray_distances.begin() + half_ray_count,
                   ray_distances.end());
  if( accepted_ray_count % 2 == 0) {
    double median_1 = ray_distances[half_ray_count];
    double median_2 = *std::max_element(ray_distances.begin(),
                                        ray_distances.begin() + half_ray_count);
    median_sdf = (median_1 + median_2) / 2;
  } else {
    median_sdf = ray_distances[half_ray_count];
  }
  /* Calculate st dev using mean_sdf as mean */
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it) {
    double dif = (*dist_it) - mean_sdf;
    st_dev += dif * dif;
  }
  st_dev = CGAL::sqrt(st_dev / (ray_distances.size()));
  /* Calculate sdf, accept rays : ray_dist - median < st dev */
  w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end();
      ++dist_it, ++w_it) {
    if(fabs((*dist_it) - median_sdf) > st_dev) {
      continue;
    }
    total_distance += (*dist_it) * (*w_it);
    total_weights += (*w_it);
  }
  return total_distance / total_weights;
}

template <class Polyhedron>
inline double
Surface_mesh_segmentation<Polyhedron>::calculate_sdf_value_from_rays_with_mean(
  std::vector<double>& ray_distances,
  std::vector<double>& ray_weights) const
{
  double total_weights = 0.0, total_distance = 0.0;
  double mean_sdf = 0.0, st_dev = 0.0;
  std::vector<double>::iterator w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it, ++w_it) {
    total_distance += (*dist_it) * (*w_it);
    total_weights += (*w_it);
  }
  mean_sdf = total_distance / total_weights;
  total_weights = 0.0;
  total_distance = 0.0;
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it) {
    double dif = (*dist_it) - mean_sdf;
    st_dev += dif * dif;
  }
  st_dev = CGAL::sqrt(st_dev / (ray_distances.size()));

  w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it, ++w_it) {
    if(fabs((*dist_it) - mean_sdf) > st_dev) {
      continue;
    }
    total_distance += (*dist_it) * (*w_it);
    total_weights += (*w_it);
  }
  return total_distance / total_weights;
}

template <class Polyhedron>
inline double
Surface_mesh_segmentation<Polyhedron>::calculate_sdf_value_from_rays_with_trimmed_mean(
  std::vector<double>& ray_distances,
  std::vector<double>& ray_weights) const
{
  std::vector<std::pair<double, double> > distances_with_weights;
  distances_with_weights.reserve(ray_distances.size());
  std::vector<double>::iterator w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it, ++w_it) {
    distances_with_weights.push_back(std::pair<double, double>((*dist_it),
                                     (*w_it)));
  }
  std::sort(distances_with_weights.begin(), distances_with_weights.end(),
            compare_pairs_using_first<std::pair<double, double> >());
  int b = floor(distances_with_weights.size() / 20.0 + 0.5); // Eliminate %5.
  int e = distances_with_weights.size() - b;                 // Eliminate %5.

  double total_weights = 0.0, total_distance = 0.0;
  double trimmed_mean = 0.0, st_dev = 0.0;

  for(int i = b; i < e; ++i) {
    total_distance += distances_with_weights[i].first *
                      distances_with_weights[i].second;
    total_weights += distances_with_weights[i].second;
  }
  trimmed_mean = total_distance / total_weights;

  total_weights = 0.0;
  total_distance = 0.0;
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it) {
    double dif = (*dist_it) - trimmed_mean;
    st_dev += dif * dif;
  }
  st_dev = CGAL::sqrt(st_dev / (ray_distances.size()));

  w_it = ray_weights.begin();
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it, ++w_it) {
    if(fabs((*dist_it) - trimmed_mean) > st_dev) {
      continue;
    }
    total_distance += (*dist_it) * (*w_it);
    total_weights += (*w_it);
  }
  return total_distance / total_weights;
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::disk_sampling_rejection()
{
  int number_of_points_sqrt = number_of_rays_sqrt;
  double length_of_normal = 1.0 / tan(cone_angle / 2);
  double mid_point = (number_of_points_sqrt-1) / 2.0;
  double angle_st_dev = cone_angle / ANGLE_ST_DEV_DIVIDER;

  for(int i = 0; i < number_of_points_sqrt; ++i)
    for(int j = 0; j < number_of_points_sqrt; ++j) {
      double w1 = (i - mid_point)/(mid_point);
      double w2 = (j - mid_point)/(mid_point);
      double R = sqrt(w1*w1 + w2*w2);
      if(R > 1.0) {
        continue;
      }
      double angle = atan(R / length_of_normal);
      angle =  exp(-0.5 * (pow(angle / angle_st_dev, 2))); // weight
      disk_samples.push_back(Disk_sample(w1, w2, angle));
    }
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::disk_sampling_polar_mapping()
{
  int number_of_points_sqrt = number_of_rays_sqrt;
  double length_of_normal = 1.0 / tan(cone_angle / 2);
  double angle_st_dev = cone_angle / ANGLE_ST_DEV_DIVIDER;

  for(int i = 0; i < number_of_points_sqrt; ++i)
    for(int j = 0; j < number_of_points_sqrt; ++j) {
      double w1 = i / (double) (number_of_points_sqrt-1);
      double w2 = j / (double) (number_of_points_sqrt-1);
      double R = w1;
      double Q = 2 * w1 * CGAL_PI;
      double angle = atan(R / length_of_normal);
      angle =  exp(-0.5 * (pow(angle / angle_st_dev, 2))); // weight
      disk_samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), angle));
    }
}

template <class Polyhedron>
inline void
Surface_mesh_segmentation<Polyhedron>::disk_sampling_concentric_mapping()
{
  int number_of_points_sqrt = number_of_rays_sqrt;
  double length_of_normal = 1.0 / tan(cone_angle / 2);
  double fraction = 2.0 / (number_of_points_sqrt -1);
  double angle_st_dev = cone_angle / ANGLE_ST_DEV_DIVIDER;

  for(int i = 0; i < number_of_points_sqrt; ++i)
    for(int j = 0; j < number_of_points_sqrt; ++j) {
      double w1 = -1 + i * fraction;
      double w2 = -1 + j * fraction;
      double R, Q;
      if(w1 == 0 && w2 == 0) {
        R = 0;
        Q = 0;
      } else if(w1 > -w2) {
        if(w1 > w2) {
          R = w1;
          Q = (w2 / w1);
        } else        {
          R = w2;
          Q = (2 - w1 / w2);
        }
      } else {
        if(w1 < w2) {
          R = -w1;
          Q = (4 + w2 / w1);
        } else        {
          R = -w2;
          Q = (6 - w1 / w2);
        }
      }
      Q *= (0.25 * CGAL_PI);
      double angle = atan(R / length_of_normal);
      angle =  exp(-0.5 * (pow(angle / angle_st_dev, 2))); // weight
      disk_samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), angle));
    }
}


template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::normalize_sdf_values()
{
  //SL: use CGAL::min_max_element //IOY: done.
  typedef typename Face_value_map::iterator fv_iterator;
  compare_pairs<typename Face_value_map::value_type> comparator;
  std::pair<fv_iterator, fv_iterator> min_max_pair =
    CGAL::min_max_element(sdf_values.begin(), sdf_values.end(), comparator,
                          comparator);

  double max_value = min_max_pair.second->second,
         min_value = min_max_pair.first->second;
  double max_min_dif = max_value - min_value;
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    double linear_normalized = (pair_it->second - min_value) / max_min_dif;
    double log_normalized = log(linear_normalized * NORMALIZATION_ALPHA + 1) /
                            LOG_5;
    pair_it->second = log_normalized;
  }
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::smooth_sdf_values()
{
  Face_value_map smoothed_sdf_values;
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    Facet_handle f = pair_it->first;
    typename Facet::Halfedge_around_facet_circulator facet_circulator =
      f->facet_begin();
    double total_neighbor_sdf = 0.0;
    do {
      total_neighbor_sdf += sdf_values[facet_circulator->opposite()->facet()];
    } while( ++facet_circulator !=  f->facet_begin());
    total_neighbor_sdf /= 3;
    smoothed_sdf_values[f] = (sdf_values[f] + total_neighbor_sdf) / 2;
  }
  sdf_values = smoothed_sdf_values;
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::apply_GMM_fitting()
{
  centers.clear();
  std::vector<double> sdf_vector;
  sdf_vector.reserve(sdf_values.size());
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    sdf_vector.push_back(pair_it->second);
  }
  Expectation_maximization fitter(number_of_centers, sdf_vector);
  std::vector<int> center_memberships;
  fitter.fill_with_center_ids(center_memberships);
  std::vector<int>::iterator center_it = center_memberships.begin();
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it, ++center_it) {
    centers.insert(std::pair<Facet_handle, int>(pair_it->first, (*center_it)));
  }
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::apply_K_means_clustering()
{
  centers.clear();
  std::vector<double> sdf_vector;
  sdf_vector.reserve(sdf_values.size());
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    sdf_vector.push_back(pair_it->second);
  }
  K_means_clustering clusterer(number_of_centers, sdf_vector);
  std::vector<int> center_memberships;
  clusterer.fill_with_center_ids(center_memberships);
  std::vector<int>::iterator center_it = center_memberships.begin();
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it, ++center_it) {
    centers.insert(std::pair<Facet_handle, int>(pair_it->first, (*center_it)));
  }
  //center_memberships_temp = center_memberships; //remove
}
template <class Polyhedron>
inline void
Surface_mesh_segmentation<Polyhedron>::apply_GMM_fitting_with_K_means_init()
{
  centers.clear();
  std::vector<double> sdf_vector;
  sdf_vector.reserve(sdf_values.size());
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    sdf_vector.push_back(pair_it->second);
  }
  K_means_clustering clusterer(number_of_centers, sdf_vector);
  std::vector<int> center_memberships;
  clusterer.fill_with_center_ids(center_memberships);
  //std::vector<int> center_memberships = center_memberships_temp;
  Expectation_maximization fitter(number_of_centers, sdf_vector,
                                  center_memberships);
  center_memberships.clear();
  fitter.fill_with_center_ids(center_memberships);
  std::vector<int>::iterator center_it = center_memberships.begin();
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it, ++center_it) {
    centers.insert(std::pair<Facet_handle, int>(pair_it->first, (*center_it)));
  }
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::write_sdf_values(
  const char* file_name)
{
  std::ofstream output(file_name);
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    output << sdf_values[facet_it] << std::endl;
  }
  output.close();
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::read_sdf_values(
  const char* file_name)
{
  std::ifstream input(file_name);
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    double sdf_value;
    input >> sdf_value;
    sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf_value));
  }
}
} //namespace CGAL
#undef ANGLE_ST_DEV_DIVIDER
#undef LOG_5
#undef NORMALIZATION_ALPHA
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H