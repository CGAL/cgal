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
private:
	// member data
	Polyhedron *m_pPolyhedron;
	std::list<Point> m_points;
	std::list<Segment> m_segments;

	// distance functions (simple 2D arrays)
	Color_ramp m_red_ramp;
	Color_ramp m_blue_ramp;
	Color_ramp m_thermal_ramp;
	FT m_max_distance_function;
	bool m_view_distance_function;
	bool m_signed_distance_function;
	typedef std::pair<Point,FT> Point_distance;
	Point_distance m_distance_function[100][100];

public:
	Scene();
	~Scene();

public:
	void draw();

	typedef CGAL::Bbox_3 Bbox;
	Bbox bbox() { return Bbox(); }

public:

	// utility functions
	Ray random_ray();
	Line random_line();
	Point random_point();
	Plane random_plane();
	Vector random_vector();
	Segment random_segment();

	// file menu
	int open(QString filename);

	// edit menu
	void clear_points() { m_points.clear(); }
	void clear_segments() { m_segments.clear(); }
	void clear_distance_function() { m_max_distance_function = 0.0; }

	// algorithms
	void generate_edge_points(const unsigned int nb_points);
	void generate_inside_points(const unsigned int nb_points);
	void generate_boundary_points(const unsigned int nb_points);
	void generate_boundary_segments(const unsigned int nb_slices);
	void refine_bisection(const FT max_sqlen);

	// distance functions 
	void signed_distance_function();
	void unsigned_distance_function();
	void unsigned_distance_function_to_edges();

	// toggle view options
	void toggle_view_points();
	void toggle_view_segments();
	void toggle_view_poyhedron();
	void toggle_view_distance_function();

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
	enum {DO_INTERSECT,
		  ANY_INTERSECTION,
		  NB_INTERSECTIONS,
		  ALL_INTERSECTIONS,
		  ALL_INTERSECTED_PRIMITIVES};
	void bench_memory();
	unsigned int nb_digits(const unsigned int value);
	void benchmark_intersections(const int duration);
	void bench_intersection_rays(Facet_tree& tree,const int function,const int duration);
	void bench_intersection_lines(Facet_tree& tree,const int function,const int duration);
	void bench_intersection_planes(Facet_tree& tree,const int function,const int duration);
	void bench_intersection_segments(Facet_tree& tree,const int function,const int duration);

	// distance benchmarks
	enum {SQ_DISTANCE,
	      CLOSEST_POINT,
	      CLOSEST_POINT_AND_PRIMITIVE_ID};
	void benchmark_distances(const int duration);
	void bench_closest_point(Facet_tree& tree,const int duration);
	void bench_squared_distance(Facet_tree& tree,const int duration);
	void bench_closest_point_and_primitive(Facet_tree& tree,const int duration);
	void bench_distance(Facet_tree& tree,const int function,const int duration);

	// intersection benchmarks
	void bench_do_intersect(Facet_tree& tree,const int duration);
	void bench_nb_intersections(Facet_tree& tree,const int duration);
	void bench_any_intersection(Facet_tree& tree,const int duration);
	void bench_all_intersections(Facet_tree& tree,const int duration);
	void bench_all_intersected_primitives(Facet_tree& tree,const int duration);

	// drawing
	void draw_points();
	void draw_segments();
	void draw_polyhedron();
	void draw_signed_distance_function();
	void draw_unsigned_distance_function();
}; // end class Scene


#endif // SCENE_H
