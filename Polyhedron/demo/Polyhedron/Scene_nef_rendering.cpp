#include "Nef_type.h"
#include <CGAL/Nef_S2/OGL_base_object.h>
#include <CGAL/Nef_3/OGL_helper.h>

void gl_render_nef_facets(Nef_polyhedron *p)
{
  CGAL::OGL::OGL_base_object* object_ = new CGAL::OGL::Polyhedron();
  CGAL::OGL::Nef3_Converter<Nef_polyhedron>::convert_to_OGLPolyhedron(*p,
                                 static_cast<CGAL::OGL::Polyhedron*>(object_));
  object_->draw();
  delete object_;
}
