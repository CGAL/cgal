#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>

#include "types.h"

#include <iostream>
#include <cmath>

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

  // TO DEFINE
  Bbox bbox() { return Bbox(); }

public:
  int open(QString filename);

  void generate_inside_points(const unsigned int nb_trials);
  void generate_boundary_segments(const unsigned int nb_slices);

  Point random_point();
  Vector random_vector();
  void benchmark_do_intersect();

  // toggle view options
  void toggle_view_poyhedron();

private:
  Polyhedron *m_pPolyhedron;
  std::list<Point> m_points;
  std::list<Segment> m_segments;
  
  // TODO: unsigned and signed distance function (shown in cut plane)

  // view options
  bool m_view_polyhedron;

}; // end class Scene


#endif // SCENE_H
