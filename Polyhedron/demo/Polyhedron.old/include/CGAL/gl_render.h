#ifndef _GL_RENDER_
#define _GL_RENDER_

#include <CGAL/gl.h>
#include <CGAL/compute_normal.h>

template <class Polyhedron>
void gl_render_facets(Polyhedron& polyhedron)
{
  typedef typename Polyhedron::Traits	    Kernel;
  typedef typename Kernel::Point_3	    Point;
  typedef typename Kernel::Vector_3	    Vector;
  typedef typename Polyhedron::Facet	    Facet;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;

  Facet_iterator f;
  for(f = polyhedron.facets_begin();
    f != polyhedron.facets_end();
    f++)
  {
    ::glBegin(GL_POLYGON);

    // compute normal
    Vector n = compute_facet_normal<Facet,Kernel>(*f);
    ::glNormal3d(n.x(),n.y(),n.z());

    // revolve around current face to get vertices
    HF_circulator he = f->facet_begin();
    HF_circulator end = he;
    CGAL_For_all(he,end)
    {
      const Point& p = he->vertex()->point();
      ::glVertex3d(p.x(),p.y(),p.z());
    }
    ::glEnd();
  }
} // end gl_render_facets

template <class Polyhedron>
void gl_render_edges(Polyhedron& polyhedron)
{
  typedef typename Polyhedron::Traits		Kernel;
  typedef typename Kernel::Point_3		Point;
  typedef typename Polyhedron::Edge_iterator	Edge_iterator;

  ::glBegin(GL_LINES);
  Edge_iterator he;
  for(he = polyhedron.edges_begin();
    he != polyhedron.edges_end();
    he++)
  {
    const Point& a = he->vertex()->point();
    const Point& b = he->opposite()->vertex()->point();
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
  }
  ::glEnd();
} // end gl_render_edges


#endif // _GL_RENDER_




