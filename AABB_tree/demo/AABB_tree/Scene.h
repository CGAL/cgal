#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>
#include "types.h"

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

  struct Bbox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    Bbox(const double _xmin,const double _ymin,const double _zmin,
         const double _xmax,const double _ymax,const double _zmax)
	 : xmin(_xmin), ymin(_ymin), zmin(_zmin),
	   xmax(_xmax), ymax(_ymax), zmax(_zmax)
    {
    }
    Bbox()
	 : xmin(0.0), ymin(0.0), zmin(0.0),
	   xmax(1.0), ymax(1.0), zmax(1.0)
    {
    }
  };

  double len_diagonal()
  {
    Bbox box = bbox();
    double dx = box.xmax - box.xmin;
    double dy = box.ymax - box.ymin;
    double dz = box.zmax - box.zmin;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  // TODO
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

  // toggle view options
  void toggle_view_points();
  void toggle_view_segments();
  void toggle_view_poyhedron();

private:
	// member data
  Polyhedron *m_pPolyhedron;
  std::list<Point> m_points;
  std::list<Segment> m_segments;
  
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

  void bench_do_intersect(Facet_tree& tree);
  void bench_nb_intersections(Facet_tree& tree);
  void bench_any_intersection(Facet_tree& tree);
  void bench_all_intersections(Facet_tree& tree);
  void bench_all_intersected_primitives(Facet_tree& tree);
}; // end class Scene


#endif // SCENE_H
