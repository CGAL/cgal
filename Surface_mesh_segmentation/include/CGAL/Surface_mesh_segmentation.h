#ifndef CGAL_SURFACE_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#define PI 3.14159265359
#define LOG_5 1.60943791
#define NORMALIZATION_ALPHA 4.0

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
  typedef typename Kernel::FT         FT;
  //Should we access types (inner classes etc.) via traits ? may be something like this:
  //typedef typename boost::graph_traits<Polyhedron>::facet_iterator facet_iterator; ???
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle   Facet_handle;
protected:

  typedef typename Kernel::Ray_3 Ray;

  typedef typename CGAL::AABB_polyhedron_triangle_primitive<Kernel, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<CGAL::AABB_traits<Kernel, Primitive>>
      Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;

  typedef std::map<Facet_handle, FT> Face_value_map;
  template <typename ValueTypeName>
  struct compare_pairs {
    bool operator()(ValueTypeName& v1, ValueTypeName& v2) {
      return v1.second < v2.second;
    }
  };
//member variables
public:
  Polyhedron* mesh;
  Face_value_map sdf_values;
protected:
  std::ofstream log_file;

//member functions
public:
  Surface_mesh_segmentation(Polyhedron* mesh) : mesh(mesh),
    log_file("log_file.txt") {
    calculate_sdf_values();
  }
  ~Surface_mesh_segmentation() {
  }

//protected:
  void calculate_sdf_values() {
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
      FT sdf = calculate_sdf_value_of_facet(facet_it, center, normal, tree,
                                            (1.0/3.0) * PI, 6);
      sdf_values.insert(std::pair<Facet_handle, FT>(facet_it, sdf));
    }
    normalize_sdf_values();
    smooth_sdf_values();
  }

  FT calculate_sdf_value_of_facet(const Facet_handle& facet, const Point& center,
                                  const Vector& normal_const, const Tree& tree, double half_cone_angle,
                                  int ray_count_sqrt) const {
    Kernel::Plane_3 plane(center, normal_const);
    Vector v1 = plane.base1();
    Vector v2 = plane.base2();
    v1 = v1 / CGAL::sqrt(v1.squared_length());
    v2 = v2 / CGAL::sqrt(v2.squared_length());
    //Vector v1 = on_plane - center;
    //v1 = v1 / sqrt(v1.squared_length());
    //Vector v2 = CGAL::cross_product(normal, v1);

    int ray_count = ray_count_sqrt * ray_count_sqrt;
    std::vector<FT> ray_distances(ray_count), ray_weights(ray_count);
    FT total_weights = FT(0.0), total_distance = FT(0.0);
    double angle_st_dev = half_cone_angle / 2; //Not sure what to use here.
    double normal_distance = 1.0 / tan(half_cone_angle);
    Vector normal = normal_const * normal_distance;
    double mid_point = (ray_count_sqrt-1) / 2.0;
    for(int i = 0; i < ray_count_sqrt; ++i)
      for(int j = 0; j < ray_count_sqrt; ++j) {
        double picking_1 = i / (double) (ray_count_sqrt-1);
        double picking_2 = j / (double) (ray_count_sqrt-1);
        double R = picking_1;
        double Q = 2 * picking_2 * PI;
        Vector random_vector = (v1 * (R * cos(Q))) + (v2 * (R * sin(Q)));
        double dist_to_center = R;
        //double w1 = (i - mid_point)/(mid_point);
        //double w2 = (j - mid_point)/(mid_point);
        //double dist_to_center = sqrt(w1*w1 + w2*w2);
        //if(dist_to_center > 1.0) { continue; }
        //Vector random_vector = (v1 * w1) + (v2 * w2);
        //random_vector = random_vector * tan(half_cone_angle);

        Ray ray(center, normal + random_vector);
        FT min_distance;
        bool is_found;
        cast_and_return_minimum(ray, tree, facet, is_found, min_distance);
        if(!is_found) {
          continue;
        }

        double angle = atan(dist_to_center / normal_distance);
        FT weight = FT(exp(-0.5 * (pow(angle / angle_st_dev, 2))));

        ray_weights.push_back(weight);
        ray_distances.push_back(min_distance);
        total_weights  += weight;
        total_distance += (min_distance * weight);
      }

    FT average_sdf = total_distance / total_weights;
    total_weights = total_distance = FT(0.0);
    FT st_dev = FT(0.0);

    for(std::vector<FT>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end(); ++dist_it) {
      FT dif = (*dist_it) - average_sdf;
      st_dev += dif * dif;
    }
    st_dev = CGAL::sqrt(st_dev / (ray_distances.size()));

    std::vector<FT>::iterator w_it = ray_weights.begin();
    for(std::vector<FT>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end(); ++dist_it, ++w_it) {
      if(abs((*dist_it) - average_sdf) > st_dev) {
        continue;
      }
      total_distance += (*dist_it) * (*w_it);
      total_weights += (*w_it);
    }
    return total_distance / total_weights;
  }

  void cast_and_return_minimum(const Ray& ray, const Tree& tree,
                               const Facet_handle& facet,
                               bool& is_found, FT& min_distance) const {
    std::list<Object_and_primitive_id> intersections;
    tree.all_intersections(ray, std::back_inserter(intersections));
    Vector min_i_ray;
    Tree::Primitive_id min_id;
    is_found = false;
    for(std::list<Object_and_primitive_id>::iterator op_it = intersections.begin();
        op_it != intersections.end() ; ++op_it) {
      CGAL::Object object   = op_it->first;
      Tree::Primitive_id id = op_it->second;
      Point i_point;
      if(id == facet)                    {
        continue;  //Since center is located on related facet, we should skip it if there is an intersection with it.
      }
      if(!CGAL::assign(i_point, object)) {
        continue;  //What to do here (in case of intersection object is a segment), I am not sure ???
      }
      Vector i_ray = (ray.source() - i_point);
      FT new_distance = CGAL::sqrt(i_ray.squared_length());
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
      //log_file << "missing" << std::endl;
    }
  }

  void normalize_sdf_values() {
    FT max_value = std::max_element(sdf_values.begin(), sdf_values.end(),
                                    compare_pairs<Face_value_map::value_type>())->second;
    FT min_value = std::min_element(sdf_values.begin(), sdf_values.end(),
                                    compare_pairs<Face_value_map::value_type>())->second;
    FT max_min_dif = max_value - min_value;
    for(Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      FT linear_normalized = (pair_it->second - min_value) / max_min_dif;
      double log_normalized = log(CGAL::to_double(linear_normalized) *
                                  NORMALIZATION_ALPHA + 1) / LOG_5; // how to log a generic Number_type(FT)?.
      pair_it->second = FT(log_normalized);
    }
  }

  void smooth_sdf_values() {
    Face_value_map smoothed_sdf_values;
    for(Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      Facet_handle f = pair_it->first;
      Facet::Halfedge_around_facet_circulator facet_circulator = f->facet_begin();
      FT total_neighbor_sdf = FT(0.0);
      do {
        total_neighbor_sdf += sdf_values[facet_circulator->opposite()->facet()];
      } while( ++facet_circulator !=  f->facet_begin());
      total_neighbor_sdf /= 3;
      smoothed_sdf_values[f] = (sdf_values[f] + total_neighbor_sdf) / 2;
    }
    sdf_values = smoothed_sdf_values;
  }

};
} //namespace CGAL
#undef LOG_5
#undef PI
#undef NORMALIZATION_ALPHA
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H