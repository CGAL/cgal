#ifndef CGAL_SURFACE_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_H
/* NEED TO BE DONE
 * About implementation:
 * +) I am not using BGL, as far as I checked there is a progress on BGL redesign
 *    (https://cgal.geometryfactory.com/CGAL/Members/wiki/Features/BGL) which introduces some features
 *    for face-based traversal / manipulation by FaceGraphs
 * +) Deciding on which parameters will be taken from user
 *
 * About paper (and correctness / efficiency etc.):
 * +) Weighting ray distances with inverse of their angles: not sure how to weight exactly
 * +) Anisotropic smoothing: have no idea what it is exactly, should read some material (google search is not enough)
 * +) Deciding how to generate rays in cone: for now using "polar angle" and "accept-reject (square)" and "concentric mapping" techniques
 */


#include "Expectation_maximization.h"
#include "K_means_clustering.h"
#include "Timer.h"

//#include <CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h>
//#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/utility.h>

#include <boost/optional.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

#define LOG_5 1.60943791
#define NORMALIZATION_ALPHA 4.0
#define ANGLE_ST_DEV_DIVIDER 3.0

//#define ACTIVATE_SEGMENTATION_DEBUG
#ifdef ACTIVATE_SEGMENTATION_DEBUG
#define SEG_DEBUG(x) x;
#else
#define SEG_DEBUG(x)
#endif

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
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Edge_iterator   Edge_iterator;
  typedef typename Polyhedron::Vertex_handle   Vertex_handle;
protected:
  typedef typename Kernel::Ray_3   Ray;
  typedef typename Kernel::Plane_3 Plane;

  typedef typename CGAL::AABB_polyhedron_triangle_primitive<Kernel, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<CGAL::AABB_traits<Kernel, Primitive> >
  Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;
  typedef typename Tree::Primitive_id
  Primitive_id;

  typedef std::map<Facet_handle, double> Face_value_map;
  typedef std::map<Facet_handle, int>    Face_center_map;
  typedef std::map<Halfedge_handle, double> Edge_angle_map;
  /*Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = angle with cone-normal (weight). */
  typedef CGAL::Triple<double, double, double>               Disk_sample;
  typedef std::vector<CGAL::Triple<double, double, double> > Disk_samples_list;

  template <typename ValueTypeName>
  struct Compare_second_element {
    bool operator()(const ValueTypeName& v1,const ValueTypeName& v2) const {
      return v1.second < v2.second;
    }
  };

//member variables
public:
  Polyhedron* mesh;

  Face_value_map  sdf_values;
  Face_center_map centers;
  Edge_angle_map  dihedral_angles;

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
                            int number_of_rays_sqrt = 7, double cone_angle = (2.0 / 3.0) * CGAL_PI,
                            int number_of_centers = 2);

  void calculate_sdf_values();

  double calculate_sdf_value_of_facet (const Facet_handle& facet,
                                       const Tree& tree) const;
  boost::optional<double> cast_and_return_minimum(const Ray& ray,
      const Tree& tree, const Facet_handle& facet) const;
  boost::optional<double> cast_and_return_minimum_use_closest(const Ray& ray,
      const Tree& tree, const Facet_handle& facet) const;

  double calculate_sdf_value_from_rays (std::vector<double>& ray_distances,
                                        std::vector<double>& ray_weights) const;
  double calculate_sdf_value_from_rays_with_mean (std::vector<double>&
      ray_distances, std::vector<double>& ray_weights) const;
  double calculate_sdf_value_from_rays_with_trimmed_mean (
    std::vector<double>& ray_distances, std::vector<double>& ray_weights) const;

  void arrange_center_orientation(const Plane& plane, const Vector& unit_normal,
                                  Point& center) const;
  void calculate_dihedral_angles();
  double calculate_dihedral_angle_of_edge(const Halfedge_handle& edge) const;

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
  void read_center_ids(const char* file_name);
};

template <class Polyhedron>
inline Surface_mesh_segmentation<Polyhedron>::Surface_mesh_segmentation(
  Polyhedron* mesh, int number_of_rays_sqrt, double cone_angle,
  int number_of_centers)
  : mesh(mesh), cone_angle(cone_angle), number_of_rays_sqrt(number_of_rays_sqrt),
    number_of_centers(number_of_centers), log_file("log_file.txt")
{
  SEG_DEBUG(Timer t)

  disk_sampling_concentric_mapping();
  //calculate_sdf_values();
  SEG_DEBUG(std::cout << t)
  //apply_GMM_fitting();
  //write_sdf_values("sdf_values_sample_dino_ws.txt");
  //read_sdf_values("sdf_values_sample_camel.txt");
  calculate_dihedral_angles();
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::calculate_sdf_values()
{
  sdf_values.clear();
  Tree tree(mesh->facets_begin(), mesh->facets_end());
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    CGAL_precondition(facet_it->is_triangle()); //Mesh should contain triangles.

    //SL: cone angle and number of rays should be parameters.
    double sdf = calculate_sdf_value_of_facet(facet_it, tree);
    sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf));
  }
  normalize_sdf_values();
  smooth_sdf_values();
}

template <class Polyhedron>
inline double
Surface_mesh_segmentation<Polyhedron>::calculate_sdf_value_of_facet(
  const Facet_handle& facet, const Tree& tree) const
{
  Point p1 = facet->halfedge()->vertex()->point();
  Point p2 = facet->halfedge()->next()->vertex()->point();
  Point p3 = facet->halfedge()->next()->next()->vertex()->point();
  Point center  = CGAL::centroid(p1, p2, p3);
  Vector normal = CGAL::unit_normal(p1, p2,
                                    p3) * -1.0; //Assuming triangles are CCW oriented.

  Plane plane(center, normal);
  Vector v1 = plane.base1();
  Vector v2 = plane.base2();
  v1 = v1 / sqrt(v1.squared_length());
  v2 = v2 / sqrt(v2.squared_length());

  arrange_center_orientation(plane, normal, center);

  int ray_count = number_of_rays_sqrt * number_of_rays_sqrt;

  std::vector<double> ray_distances, ray_weights;
  ray_distances.reserve(ray_count);
  ray_weights.reserve(ray_count);

  double length_of_normal = 1.0 / tan(cone_angle / 2);
  normal = normal * length_of_normal;

  for(Disk_samples_list::const_iterator sample_it = disk_samples.begin();
      sample_it != disk_samples.end(); ++sample_it) {
    Vector disk_vector = v1 * sample_it->first + v2 * sample_it->second;
    Ray ray(center, normal + disk_vector);

    //boost::optional<double> min_distance = cast_and_return_minimum(ray, tree, facet);
    boost::optional<double> min_distance = cast_and_return_minimum_use_closest(ray,
                                           tree, facet);
    if(!min_distance) {
      continue;
    }

    ray_weights.push_back(sample_it->third);
    ray_distances.push_back(*min_distance);
  }
  return calculate_sdf_value_from_rays_with_trimmed_mean(ray_distances,
         ray_weights);
}

template <class Polyhedron>
boost::optional<double>
Surface_mesh_segmentation<Polyhedron>::cast_and_return_minimum(
  const Ray& ray, const Tree& tree, const Facet_handle& facet) const
{
  boost::optional<double> min_distance;
  std::list<Object_and_primitive_id> intersections;
  tree.all_intersections(ray, std::back_inserter(intersections));
  Vector min_i_ray;
  Primitive_id min_id;
  for(typename std::list<Object_and_primitive_id>::iterator op_it =
        intersections.begin();
      op_it != intersections.end() ; ++op_it) {
    CGAL::Object object = op_it->first;
    Primitive_id id     = op_it->second;
    if(id == facet) {
      continue;  //Since center is located on related facet, we should skip it if there is an intersection with it.
    }

    Point i_point;
    if(!CGAL::assign(i_point, object)) {
      continue;  //What to do here (in case of intersection object is a segment), I am not sure ???
    }

    Vector i_ray = (ray.source() - i_point);
    double new_distance = i_ray.squared_length();
    if(!min_distance || new_distance < min_distance) {
      min_distance = new_distance;
      min_id = id;
      min_i_ray = i_ray;
    }
  }
  if(!min_distance) {
    return min_distance;
  }

  Point min_v1 = min_id->halfedge()->vertex()->point();
  Point min_v2 = min_id->halfedge()->next()->vertex()->point();
  Point min_v3 = min_id->halfedge()->next()->next()->vertex()->point();
  Vector min_normal = CGAL::normal(min_v1, min_v2, min_v3) * -1.0;

  if(CGAL::angle(CGAL::ORIGIN + min_i_ray, Point(CGAL::ORIGIN),
                 CGAL::ORIGIN + min_normal) != CGAL::ACUTE) {
    return boost::optional<double>();
  }
  return (min_distance = sqrt(*min_distance));
}

template <class Polyhedron>
boost::optional<double>
Surface_mesh_segmentation<Polyhedron>::cast_and_return_minimum_use_closest (
  const Ray& ray, const Tree& tree,
  const Facet_handle& facet) const
{
  //static double dist = 0.1;
  //boost::optional<double> min_distance_2 = dist++;
  //return min_distance_2;

  boost::optional<double> min_distance;
  boost::optional<Object_and_primitive_id> intersection =
    tree.closest_intersection(ray);
  //boost::optional<Object_and_primitive_id> intersection = tree.any_intersection(ray);
  if(!intersection) {
    return min_distance;
  }

  CGAL::Object object = intersection->first;
  Primitive_id min_id = intersection->second;

  if(min_id == facet) {
    CGAL_error();  // There should be no such case, after center-facet arrangments.
  }

  Point i_point;
  if(!CGAL::assign(i_point, object)) {
    return min_distance;
  }

  Vector min_i_ray = ray.source() - i_point;
  Point min_v1 = min_id->halfedge()->vertex()->point();
  Point min_v2 = min_id->halfedge()->next()->vertex()->point();
  Point min_v3 = min_id->halfedge()->next()->next()->vertex()->point();
  Vector min_normal = CGAL::normal(min_v1, min_v2, min_v3) * -1.0;

  if(CGAL::angle(CGAL::ORIGIN + min_i_ray, Point(CGAL::ORIGIN),
                 CGAL::ORIGIN + min_normal) != CGAL::ACUTE) {
    return min_distance;
  }
  return (min_distance = sqrt(min_i_ray.squared_length()));
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
  st_dev = sqrt(st_dev / ray_distances.size());
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
  st_dev = sqrt(st_dev / ray_distances.size());

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
    distances_with_weights.push_back(std::pair<double, double>(*w_it, *dist_it));
  }
  std::sort(distances_with_weights.begin(), distances_with_weights.end(),
            Compare_second_element<std::pair<double, double> >());
  int b = floor(distances_with_weights.size() / 20.0 + 0.5); // Eliminate %5.
  int e = distances_with_weights.size() - b;                 // Eliminate %5.

  double total_weights = 0.0, total_distance = 0.0;
  double trimmed_mean = 0.0, st_dev = 0.0;

  for(int i = b; i < e; ++i) {
    total_distance += distances_with_weights[i].first *
                      distances_with_weights[i].second;
    total_weights += distances_with_weights[i].first;
  }
  trimmed_mean = total_distance / total_weights;

  total_weights = 0.0;
  total_distance = 0.0;
  for(std::vector<double>::iterator dist_it = ray_distances.begin();
      dist_it != ray_distances.end(); ++dist_it) {
    double dif = (*dist_it) - trimmed_mean;
    st_dev += dif * dif;
  }
  st_dev = sqrt(st_dev / ray_distances.size());

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

/* Slightly move center towards inverse normal of facet.
 * Parameter plane is constructed by inverse normal.
 */
template <class Polyhedron>
void Surface_mesh_segmentation<Polyhedron>::arrange_center_orientation(
  const Plane& plane, const Vector& unit_normal, Point& center) const
{
  /*
  *  Not sure how to decide on how much I should move center ?
  */
  double epsilon = 1e-8;
  // Option-1
  //if(plane.has_on_positive_side(center)) { return; }
  //Vector epsilon_normal = unit_normal * epsilon;
  //do {
  //    center = CGAL::operator+(center, epsilon_normal);
  //} while(!plane.has_on_positive_side(center));

  //Option-2
  //if(plane.has_on_positive_side(center)) { return; }
  //double distance = sqrt(CGAL::squared_distance(plane, center));
  //distance = distance > epsilon ? (distance + epsilon) : epsilon;
  //Vector distance_normal = unit_normal * distance;
  //
  //do {
  //    center = CGAL::operator+(center, distance_normal);
  //} while(!plane.has_on_positive_side(center));

  //Option-3
  Ray ray(center, unit_normal);
  Kernel::Do_intersect_3 intersector = Kernel().do_intersect_3_object();
  if(!intersector(ray, plane)) {
    return;
  }

  Vector epsilon_normal = unit_normal * epsilon;
  do {
    center = CGAL::operator+(center, epsilon_normal);
    ray = Ray(center, unit_normal);
  } while(intersector(ray, plane));

}

template <class Polyhedron>
void Surface_mesh_segmentation<Polyhedron>::calculate_dihedral_angles()
{
  for(Edge_iterator edge_it = mesh->edges_begin(); edge_it != mesh->edges_end();
      ++edge_it) {
    double angle = calculate_dihedral_angle_of_edge(edge_it);
    dihedral_angles.insert(std::pair<Halfedge_handle, double>(edge_it, angle));
  }
}

// if concave then returns angle value between [epsilon - 1] which corresponds to angle [0 - Pi]
// if convex then returns epsilon directly.
template <class Polyhedron>
double Surface_mesh_segmentation<Polyhedron>::calculate_dihedral_angle_of_edge(
  const Halfedge_handle& edge) const
{
  double epsilon = 1e-8; // not sure but should not return zero for log(angle)...
  Facet_handle f1 = edge->facet();
  Facet_handle f2 = edge->opposite()->facet();

  Point f2_v1 = f2->halfedge()->vertex()->point();
  Point f2_v2 = f2->halfedge()->next()->vertex()->point();
  Point f2_v3 = f2->halfedge()->next()->next()->vertex()->point();
  /**
   * As far as I see from results, segment boundaries are occurred in 'concave valleys'.
   * There is no such thing written (clearly) in the paper but should we just penalize 'concave' edges (not convex edges) ?
   * Actually that is what I understood from 'positive dihedral angle'.
   */
  Point unshared_point_on_f1 = edge->next()->vertex()->point();
  Plane p2(f2_v1, f2_v2, f2_v3);
  bool concave = p2.has_on_positive_side(unshared_point_on_f1);
  if(!concave) {
    return epsilon;  // So no penalties for convex dihedral angle ? Not sure...
  }

  Point f1_v1 = f1->halfedge()->vertex()->point();
  Point f1_v2 = f1->halfedge()->next()->vertex()->point();
  Point f1_v3 = f1->halfedge()->next()->next()->vertex()->point();
  Vector f1_normal = CGAL::unit_normal(f1_v1, f1_v2, f1_v3);
  Vector f2_normal = CGAL::unit_normal(f2_v1, f2_v2, f2_v3);

  double dot = f1_normal * f2_normal;
  if(dot > 1.0)       {
    dot = 1.0;
  } else if(dot < -1.0) {
    dot = -1.0;
  }
  double angle = acos(dot) / CGAL_PI; // [0-1] normalize
  if(angle < epsilon) {
    angle = epsilon;
  }
  return angle;
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
  Compare_second_element<typename Face_value_map::value_type> comparator;
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
    smoothed_sdf_values[f] = (pair_it->second + total_neighbor_sdf) / 2;
  }
  sdf_values = smoothed_sdf_values;
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::apply_GMM_fitting()
{
  SEG_DEBUG(Timer tg)
  centers.clear();
  std::vector<double> sdf_vector;
  sdf_vector.reserve(sdf_values.size());
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it) {
    sdf_vector.push_back(pair_it->second);
  }
  SEG_DEBUG(Timer t)
  internal::Expectation_maximization fitter(number_of_centers, sdf_vector, 50);
  SEG_DEBUG(std::cout << t)
  std::vector<int> center_memberships;
  fitter.fill_with_center_ids(center_memberships);
  std::vector<int>::iterator center_it = center_memberships.begin();
  for(typename Face_value_map::iterator pair_it = sdf_values.begin();
      pair_it != sdf_values.end(); ++pair_it, ++center_it) {
    centers.insert(std::pair<Facet_handle, int>(pair_it->first, (*center_it)));
  }
  SEG_DEBUG(std::cout << tg)
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
  internal::K_means_clustering clusterer(number_of_centers, sdf_vector);
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
  internal::K_means_clustering clusterer(number_of_centers, sdf_vector);
  std::vector<int> center_memberships;
  clusterer.fill_with_center_ids(center_memberships);
  //std::vector<int> center_memberships = center_memberships_temp;
  internal::Expectation_maximization fitter(number_of_centers, sdf_vector,
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
  sdf_values.clear();
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    double sdf_value;
    input >> sdf_value;
    sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf_value));
  }
}

template <class Polyhedron>
inline void Surface_mesh_segmentation<Polyhedron>::read_center_ids(
  const char* file_name)
{
  std::ifstream input(file_name);
  centers.clear();
  int max_center = 0;
  for(Facet_iterator facet_it = mesh->facets_begin();
      facet_it != mesh->facets_end(); ++facet_it) {
    int center_id;
    input >> center_id;
    centers.insert(std::pair<Facet_handle, int>(facet_it, center_id));
    if(center_id > max_center) {
      max_center = center_id;
    }
  }
  number_of_centers = max_center + 1;
}
} //namespace CGAL
#undef ANGLE_ST_DEV_DIVIDER
#undef LOG_5
#undef NORMALIZATION_ALPHA

#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H