#ifndef _GL_RENDER_
#define _GL_RENDER_

#include <CGAL/gl.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>



inline void CGALglcolor(QColor c, int dv = 0)
{
  if ( 0 != dv )
  {
    c = c.darker(dv);
#undef darker
  }
  ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
}

template <class Polyhedron>
void gl_render_facets(Polyhedron& polyhedron, const std::vector<QColor>& colors)
{
  typedef typename Polyhedron::Traits	    Kernel;
  typedef typename Kernel::Point_3	    Point;
  typedef typename Kernel::Vector_3	    Vector;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;

  // Get current shading model
  GLint shading;
  ::glGetIntegerv(GL_SHADE_MODEL, &shading);

  int patch_id = -1;

  Facet_iterator f;
  for(f = polyhedron.facets_begin();
    f != polyhedron.facets_end();
    f++)
  {
    const int this_patch_id = f->patch_id();
    if(patch_id != this_patch_id) {
      CGALglcolor(colors[this_patch_id]);
      patch_id = this_patch_id;
    }
    ::glBegin(GL_POLYGON);

    // If Flat shading: 1 normal per polygon
    if (shading == GL_FLAT)
    {
      Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f, polyhedron);
      ::glNormal3d(n.x(),n.y(),n.z());
    }

    // revolve around current face to get vertices
    HF_circulator he = f->facet_begin();
    HF_circulator end = he;
    CGAL_For_all(he,end)
    {
      // If Gouraud shading: 1 normal per vertex
      if (shading == GL_SMOOTH)
      {
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(he->vertex(), polyhedron);
        ::glNormal3d(n.x(),n.y(),n.z());
      }

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




