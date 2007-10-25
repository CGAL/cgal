#ifndef _POLYHEDRAL_SURFACE_H
#define _POLYHEDRAL_SURFACE_H

#define CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

#include "surface.h"
#include "viewer.h"

// CGAL
#include <CGAL/basic.h>

// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// surface
#include <CGAL/Polyhedral_surface_3.h>

// piece-wise loop subdivision
#include <CGAL/pws_loop_subdivision.h>

#include <fstream>
#include <iostream>

#include <QAction>
#include <QMainWindow>
#include <QStatusBar>

class Polyhedral_surface : public Surface
{
  Q_OBJECT;
  
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::FT FT;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Sphere_3 Sphere;
  typedef Kernel::Vector_3 Vector;
  typedef Kernel::Triangle_3 Triangle_3;

  typedef CGAL::Polyhedral_surface_3<Kernel,
    CGAL::Surface_mesher::Has_edges> CGAL_polyhedral_surface;
  typedef CGAL_polyhedral_surface::Polyhedron Polyhedron;
public:
  Polyhedral_surface(QObject* parent,
		     double sharp_edges_angle_lower_bound,
		     double sharp_edges_angle_upper_bound)
    : Surface(parent), 
      surface_ptr(0),
      parent(parent),
      display_octree(false),
      display_surface(true),
      display_all_edges(false),
      sharp_edges_angle_lower_bound(sharp_edges_angle_lower_bound),
      sharp_edges_angle_upper_bound(sharp_edges_angle_upper_bound)
  {
    QAction* action = parent->findChild<QAction*>("actionDisplay_octree");
    if(action)
      connect(action, SIGNAL(toggled(bool)),
              this, SLOT(toggle_display_octree(bool)));

    action = parent->findChild<QAction*>("actionDisplay_surface");
    if(action)
      connect(action, SIGNAL(toggled(bool)),
              this, SLOT(toggle_display_surface(bool)));

    action = parent->findChild<QAction*>("actionDisplay_all_edges");
    if(action)
      connect(action, SIGNAL(toggled(bool)),
              this, SLOT(toggle_display_all_edges(bool)));

    action = parent->findChild<QAction*>("actionInverse_normals");
    if(action)
      connect(action, SIGNAL(toggled(bool)),
              this, SLOT(set_inverse_normals(bool)));

    action = parent->findChild<QAction*>("actionSubdivision");
    if(action)
      connect(action, SIGNAL(triggered()), this, SLOT(make_one_subdivision_step()));

    connect(parent, SIGNAL(new_sharp_edges_angle_bounds(double, double)),
            this, SLOT(set_sharp_edges_angle_bounds(double, double)));

    connect(this, SIGNAL(changed()), this, SLOT(set_message()));
  }
public slots:
  void set_message() const 
  {
    QMainWindow* mw = qobject_cast<QMainWindow *>(parent);
    if(surface_ptr && mw)
    {
      mw->statusBar()->showMessage(QString("%1 vertices. %2 edges. %3 facets.")
                                   .arg(surface_ptr->size_of_vertices())
                                   .arg(surface_ptr->size_of_halfedges()/2)
                                   .arg(surface_ptr->size_of_facets()));
    }
  }

void set_sharp_edges_angle_bounds(const double lower_bound, const double upper_bound)
  {
    sharp_edges_angle_lower_bound = lower_bound;
    sharp_edges_angle_upper_bound = upper_bound;
    surface_ptr->set_sharp_edges_angle_lower_bound(lower_bound);
    surface_ptr->set_sharp_edges_angle_upper_bound(upper_bound);
    emit changed();
  }

  void toggle_display_octree(bool b)
  {
    display_octree = b;
    emit changed();
  }

  void toggle_display_surface(bool b)
  {
    display_surface = b;
    emit changed();
  }

  void toggle_display_all_edges(bool b)
  {
    display_all_edges = b;
    emit changed();
  }

  void make_one_subdivision_step()
  {
    if(surface_ptr)
    {
      Polyhedron output;
      CSubdivider_loop<Polyhedron , Kernel> pw_loop_subdiviser;

      pw_loop_subdiviser.subdivide(*surface_ptr, output);
      static_cast<Polyhedron&>(*surface_ptr) = output;
      surface_ptr->compute_normals();
      surface_ptr->compute_type();
      emit changed();
    }
  }

public:
  void open(const QString& filename)
  {
    std::cerr << "Opening file \"" << filename.toLatin1().data() << "\"...";
    if(surface_ptr)
      delete surface_ptr;
    std::ifstream in(filename.toLatin1().data());
    surface_ptr = new CGAL_polyhedral_surface(in, 
					      sharp_edges_angle_lower_bound, 
					      sharp_edges_angle_upper_bound);
    std::cerr << " Done.\n";
    emit changed();
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
      if(display_surface)
      {
      // enable polygon offset
        ::glEnable(GL_POLYGON_OFFSET_FILL);
        ::glPolygonOffset(1.0f,1.0f);
        ::glEnable(GL_LIGHTING);
        ::glColor3f(0.2f, 0.2f, 1.f);
        surface_ptr->gl_draw_direct_triangles(false,true,inverse_normals());

        ::glDisable(GL_LIGHTING);
        ::glLineWidth(1.0f);
        if(display_all_edges)
        {
          // superimpose ordinary edges
          ::glColor3d(0.4,0.4,1.);
          surface_ptr->superimpose_edges(false,true);
        }
        // superimpose control edges
        ::glDisable(GL_LIGHTING);
        ::glColor3ub(0,0,0);
        ::glLineWidth(1.0f);
        surface_ptr->superimpose_edges(true,false);

        // draw sharp edges
        ::glColor3ub(128,128,128);
        ::glLineWidth(1.0f);
        surface_ptr->gl_draw_sharp_edges(3.0f,255,0,0);
        ::glLineWidth(1.0f);
      }
      if(display_octree)
      {
        ::glColor3ub(0,0,0);
        ::glDisable(GL_LIGHTING);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        ::glEnable(GL_LINE_STIPPLE);
        surface_ptr->gl_draw_facet_octree();
        ::glDisable(GL_LINE_STIPPLE);
        ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        ::glEnable(GL_LIGHTING);
      }

      ::glDisable(GL_POLYGON_OFFSET_FILL);
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
  QObject* parent;
  bool display_octree;
  bool display_surface;
  bool display_all_edges;
  double sharp_edges_angle_lower_bound;
  double sharp_edges_angle_upper_bound;
};

#endif // _POLYHEDRAL_SURFACE_H
