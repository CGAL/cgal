#include "Polyhedron_type.h"
#include <QtOpenGL/qgl.h>
#include <CGAL/gl_render.h>

void gl_render_facets(Polyhedron* poly)
{
  gl_render_facets(*poly);
}

void gl_render_edges(Polyhedron *poly)
{
  gl_render_edges(*poly);
}
