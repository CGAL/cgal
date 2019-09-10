#ifndef _GL_RENDER_EDGES_
#define _GL_RENDER_EDGES_

#include "Viewer.h"

template <class Polyhedron>
void gl_render_edges(Polyhedron& polyhedron, Viewer* viewer)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Kernel::Point_3 Point;

  viewer->glBegin(GL_LINES);
  typename Polyhedron::Edge_iterator he;
  for(he = polyhedron.edges_begin();
      he != polyhedron.edges_end();
      he++)
  {
    const Point& a = he->vertex()->point();
    const Point& b = he->opposite()->vertex()->point();
    viewer->glVertex3d(a.x(),a.y(),a.z());
    viewer->glVertex3d(b.x(),b.y(),b.z());
  }
  viewer->glEnd();
} 

#endif // _GL_RENDER_EDGES_




