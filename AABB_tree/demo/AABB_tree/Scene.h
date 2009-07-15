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
public:
    // types
    typedef CGAL::Bbox_3 Bbox;
    typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Facet_tree;
    typedef Facet_tree::Object_and_primitive_id Object_and_primitive_id;
    typedef Facet_tree::Primitive_id Primitive_id;

public:
    void draw(); 
    void update_bbox();
    Bbox bbox() { return m_bbox; }

private:
    // member data
    Bbox m_bbox;
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

private:
    // utility functions
    Vector random_vector();
    Ray random_ray(const Bbox& bbox);
    Line random_line(const Bbox& bbox);
    Point random_point(const Bbox& bbox);
    Plane random_plane(const Bbox& bbox);
    Segment random_segment(const Bbox& bbox);
    FT random_in(const double a,const double b);

public:
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
    void generate_points_in(const unsigned int nb_points,
        const double min, const double max);

    // algorithms/refine
    void refine_loop();
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

    // benchmarks 
    enum {DO_INTERSECT,
        ANY_INTERSECTION,
        NB_INTERSECTIONS,
        ALL_INTERSECTIONS,
        ANY_INTERSECTED_PRIMITIVE,
        ALL_INTERSECTED_PRIMITIVES};
    void bench_memory();
    void bench_construction();
    void bench_distances_vs_nbt();
    void bench_intersections_vs_nbt();
    void benchmark_intersections(const double duration);
    unsigned int nb_digits(const unsigned int value);

    template <class Query>
    void bench_intersection(Facet_tree& tree,const int function,const double duration,
        const char *query_name, const std::vector<Query>& queries, const int nb_queries);
    void bench_intersections(Facet_tree& tree, const double duration, const int function,
        const char *function_name, const std::vector<Ray>& rays,
        const std::vector<Line>& lines, const std::vector<Plane>& planes,
        const std::vector<Segment>& segments, const int nb_queries);

    // distance benchmarks
    enum {SQ_DISTANCE,
        CLOSEST_POINT,
        CLOSEST_POINT_AND_PRIMITIVE_ID};
    void benchmark_distances(const double duration);
    void bench_closest_point(Facet_tree& tree,const double duration);
    void bench_squared_distance(Facet_tree& tree,const double duration);
    void bench_closest_point_and_primitive(Facet_tree& tree,const double duration);
    void bench_distance(Facet_tree& tree,const int function,const double duration);

    // drawing
    void draw_points();
    void draw_segments();
    void draw_polyhedron();
    void draw_signed_distance_function();
    void draw_unsigned_distance_function();
}; // end class Scene


#endif // SCENE_H
