#include "Polyhedron_type.h"
#include <QtOpenGL/qgl.h>
#include <CGAL/gl_render.h>

void gl_render_polyhedron_facets(Polyhedron *p)
{
  gl_render_facets(*p);
}

void gl_render_polyhedron_edges(Polyhedron *p)
{
  gl_render_edges(*p);
}

