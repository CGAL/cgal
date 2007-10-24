#ifndef _POLYHEDRAL_SURFACE_H
#define _POLYHEDRAL_SURFACE_H

#define CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

#include "surface.h"

// CGAL
#include <CGAL/basic.h>

// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// surface
#include <CGAL/Polyhedral_surface_3.h>

#include <fstream>
#include <iostream>

class Polyhedral_surface : public Surface
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::FT FT;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Sphere_3 Sphere;
  typedef Kernel::Vector_3 Vector;
  typedef Kernel::Triangle_3 Triangle_3;

  typedef CGAL::Polyhedral_surface_3<Kernel,
    CGAL::Surface_mesher::Has_edges> CGAL_polyhedral_surface;
public:
  Polyhedral_surface() :
    surface_ptr(0) 
  {
  }

public:
  void open(const QString& filename)
  {
    std::cerr << "Opening file \"" << filename.toLatin1().data() << "\"...";
    if(surface_ptr)
      delete surface_ptr;
    std::ifstream in(filename.toLatin1().data());
    surface_ptr = new CGAL_polyhedral_surface(in, FT(0.5));
    std::cerr << " Done.\n";
  }

  void close() 
  {
    delete surface_ptr;
    surface_ptr = 0;
  }
  
  void draw() 
  {
    if(surface_ptr)
    {
      // enable polygon offset
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.0f,1.0f);
      ::glEnable(GL_LIGHTING);
      glColor3f(0.2f, 0.2f, 1.f);
      surface_ptr->gl_draw_direct_triangles(false,true);

      // superimpose edges
      ::glDisable(GL_LIGHTING);
      ::glColor3ub(0,0,0);
      ::glLineWidth(1.0f);
      surface_ptr->superimpose_edges(true,false);

      // draw sharp edges
      ::glColor3ub(128,128,128);
      ::glLineWidth(1.0f);
      surface_ptr->gl_draw_sharp_edges(3.0f,255,0,0);
      ::glLineWidth(1.0f);
    
      glDisable(GL_POLYGON_OFFSET_FILL);
    }
  }
  void get_bbox(float& xmin, float& ymin, float& zmin,
		float& xmax, float& ymax, float& zmax)
  {
    if(surface_ptr) {
      xmin=surface_ptr->bbox().xmin();
      ymin=surface_ptr->bbox().ymin();
      ymin=surface_ptr->bbox().zmin();
      xmax=surface_ptr->bbox().xmax();
      ymax=surface_ptr->bbox().ymax();
      zmax=surface_ptr->bbox().zmax();
    }
    else 
    {
      xmin = ymin = zmin = 0.f;
      xmax = ymax = zmax = 1.f;
    }
  }
    
private:
  CGAL_polyhedral_surface* surface_ptr;
};

#endif // _POLYHEDRAL_SURFACE_H
