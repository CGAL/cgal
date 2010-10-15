#include "Textured_polyhedron_type.h"
#include <QtOpenGL/qgl.h>
#include <CGAL/gl_render.h>

void gl_render_tex_polyhedron_facets(Textured_polyhedron *p)
{
  p->gl_draw_textured_triangles(true,true,1.0);
}
