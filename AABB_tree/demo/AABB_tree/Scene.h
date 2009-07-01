#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include <CGAL/AABB_intersections.h> 
#include "types.h"
#include "Color_ramp.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

class Scene
{
public:
  Scene();
  ~Scene();

  Polyhedron* polyhedron() const;

  void draw();

  typedef CGAL::Bbox_3 Bbox;
  Bbox bbox() { return Bbox(); }

public:

  // utility functions
  Point random_point();
  Vector random_vector();

  // file menu functions
  int open(QString filename);

  // edit menus functions
  void clear_points() { m_points.clear(); }
  void clear_segments() { m_segments.clear(); }

  // benchmark menu functions
  void benchmark_intersections();
  void benchmark_distances();

  // algorithms
  void generate_edge_points(const unsigned int nb_points);
  void generate_inside_points(const unsigned int nb_trials);
  void generate_boundary_points(const unsigned int nb_points);
  void generate_boundary_segments(const unsigned int nb_slices);

  // distance functions 
  void signed_distance_function();
  void unsigned_distance_function();
  void unsigned_distance_function_to_edges();

  // toggle view options
  void toggle_view_points();
  void toggle_view_segments();
  void toggle_view_poyhedron();

private:
	// member data
  Polyhedron *m_pPolyhedron;
  std::list<Point> m_points;
  std::list<Segment> m_segments;

  // distance functions (simple 2D arrays)
  FT m_max_distance_function;
  bool m_signed_distance_function;
  typedef std::pair<Point,FT> Point_distance;
  Point_distance m_distance_function[100][100];
  Color_ramp m_thermal_ramp;
  Color_ramp m_red_ramp;
  Color_ramp m_blue_ramp;
  
  // view options
  bool m_view_points;
  bool m_view_segments;
  bool m_view_polyhedron;

  // types
  typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Facet_tree;
  typedef Facet_tree::Object_and_primitive_id Object_and_primitive_id;
  typedef Facet_tree::Primitive_id Primitive_id;

  // benchmarks
  void bench_do_intersect(Facet_tree& tree);
  void bench_closest_point(Facet_tree& tree);
  void bench_squared_distance(Facet_tree& tree);
  void bench_nb_intersections(Facet_tree& tree);
  void bench_any_intersection(Facet_tree& tree);
  void bench_all_intersections(Facet_tree& tree);
  void bench_closest_point_and_primitive(Facet_tree& tree);
  void bench_all_intersected_primitives(Facet_tree& tree);

  // drawing
  void draw_points();
  void draw_segments();
  void draw_polyhedron();
  void draw_signed_distance_function();
  void draw_unsigned_distance_function();
}; // end class Scene


#endif // SCENE_H
