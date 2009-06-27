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

private:
  Polyhedron *m_pPolyhedron;
}; // end class Scene


#endif // SCENE_H
