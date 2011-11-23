#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_NEF
#include "Nef_type.h"
#include <CGAL/Nef_S2/OGL_base_object.h>
#include <CGAL/Nef_3/OGL_helper.h>
#include <CGAL/glu.h>

#include <QList>

#define GL_MACRO(type) case type: std::cerr << #type; break;

inline void CGAL_GLU_TESS_CALLBACK beginCallback(GLenum which)
{ 
//   std::cerr << "glBegin(";
//   switch(which)
//   {
//     GL_MACRO(GL_POINTS);
//     GL_MACRO(GL_LINES);
//     GL_MACRO(GL_LINE_STRIP);
//     GL_MACRO(GL_LINE_LOOP);
//     GL_MACRO(GL_TRIANGLES);
//     GL_MACRO(GL_TRIANGLE_STRIP);
//     GL_MACRO(GL_TRIANGLE_FAN);
//     GL_MACRO(GL_QUADS);
//     GL_MACRO(GL_QUAD_STRIP);
//     GL_MACRO(GL_POLYGON);
//   }
//   std::cerr << ")\n";
  glBegin(which);
}

inline void CGAL_GLU_TESS_CALLBACK endCallback(void)
{ 
//   std::cerr << "glEnd()\n";
  glEnd();
}

inline void CGAL_GLU_TESS_CALLBACK errorCallback(GLenum errorCode)
{ const GLubyte *estring;
  estring = gluErrorString(errorCode);
  fprintf(stderr, "Tessellation Error: %s\n", estring);
  std::exit (0);
}

inline void CGAL_GLU_TESS_CALLBACK vertexCallback(GLvoid* vertex,
						  GLvoid* user)
{ GLdouble* pc(static_cast<GLdouble*>(vertex));
  GLdouble* pu(static_cast<GLdouble*>(user));
  //    CGAL_NEF_TRACEN("vertexCallback coord  "<<pc[0]<<","<<pc[1]<<","<<pc[2]);
  //    CGAL_NEF_TRACEN("vertexCallback normal "<<pu[0]<<","<<pu[1]<<","<<pu[2]);
  glNormal3dv(pu);
  glVertex3dv(pc); 
//   std::cerr << "  glNormal("
// 	    << pu[0] << ", "
// 	    << pu[1] << ", "
// 	    << pu[2] << ")\n";
//   std::cerr << "  glVertex("
// 	    << pc[0] << ", "
// 	    << pc[1] << ", "
// 	    << pc[2] << ")\n";
}

struct DPoint {
  DPoint(GLdouble x, GLdouble y, GLdouble z)
  {
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
  }
  GLdouble coords[3];
};

void gl_render_nef_facets(Nef_polyhedron *p)
{
//   glPointSize(10);
//   CGAL_forall_vertices(v, p->sncp())

  
  GLboolean old_light_model;
  ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &old_light_model);

  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLUtesselator* tess_ = ::gluNewTess();
  ::gluTessCallback(tess_, GLenum(GLU_TESS_VERTEX_DATA),
		    (GLvoid (CGAL_GLU_TESS_CALLBACK *)(CGAL_GLU_TESS_DOTS)) &vertexCallback);
  ::gluTessCallback(tess_, GLenum(GLU_TESS_BEGIN),
		    (GLvoid (CGAL_GLU_TESS_CALLBACK *)(CGAL_GLU_TESS_DOTS)) &beginCallback);
  ::gluTessCallback(tess_, GLenum(GLU_TESS_END),
		    (GLvoid (CGAL_GLU_TESS_CALLBACK *)(CGAL_GLU_TESS_DOTS)) &endCallback);
  ::gluTessCallback(tess_, GLenum(GLU_TESS_ERROR),
		    (GLvoid (CGAL_GLU_TESS_CALLBACK *)(CGAL_GLU_TESS_DOTS)) &errorCallback);
  ::gluTessProperty(tess_, GLenum(GLU_TESS_WINDING_RULE),
		    GLU_TESS_WINDING_POSITIVE);

  for(Nef_polyhedron::Halffacet_const_iterator
	f = p->halffacets_begin (),
	end = p->halffacets_end();
      f != end; ++f)
  {
    if(f->is_twin()) continue;

    Nef_polyhedron::Vector_3 v = f->plane().orthogonal_vector();
    GLdouble normal[3];
    normal[0] = CGAL::to_double(v.x());
    normal[1] = CGAL::to_double(v.y());
    normal[2] = CGAL::to_double(v.z());
    GLdouble norm = normal[0]*normal[0]
      + normal[1]*normal[1]
      + normal[2]*normal[2];
    norm = CGAL::sqrt(norm);

    normal[0] /= norm;
    normal[1] /= norm;
    normal[2] /= norm;

    ::gluTessBeginPolygon(tess_, normal);
    ::gluTessNormal(tess_,normal[0],normal[1],normal[2]);

    QList<DPoint> points;
    for(Nef_polyhedron::Halffacet_cycle_const_iterator 
	  fc = f->facet_cycles_begin(),
	  end = f->facet_cycles_end();
	fc != end; ++fc)
    {
      if ( fc.is_shalfedge() )
      {
	::gluTessBeginContour(tess_);

	Nef_polyhedron::SHalfedge_const_handle h = fc;
	Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(h), he(hc);
	CGAL_For_all(hc,he){ // all vertex coordinates in facet cycle
	  Nef_polyhedron::SVertex_const_handle v = hc->source();
	  const Nef_polyhedron::Point_3& point = v->source()->point();
	  int i = points.size();
	  DPoint dp(CGAL::to_double(point.x()),
                    CGAL::to_double(point.y()),
                    CGAL::to_double(point.z()));
	  points.push_back(dp);

	  ::gluTessVertex(tess_, 
			  static_cast<GLdouble*>(static_cast<void*>(&(points[i].coords))),
			  &(points[i].coords));
	} // end facet cycles verticeses

	::gluTessEndContour(tess_);
      }
    } // end facet cycles
    ::gluTessEndPolygon(tess_);
  } // end facets
  ::gluDeleteTess(tess_);

  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, old_light_model);

  GLenum gl_error = ::glGetError();
  if(gl_error != GL_NO_ERROR)
  {
    std::cerr << "OPENGL ERROR in gl_render_nef_facets!\n  "
	      << ::gluErrorString(gl_error) << "\n";
  }
} // end gl_render_nef_facets

void gl_render_nef_edges(Nef_polyhedron *p)
{
  ::glBegin(GL_LINES);
  
  for(Nef_polyhedron::Halfedge_const_iterator 
	e = p->halfedges_begin(),
	end = p->halfedges_end();
      e != end; ++e)
  {
    const Nef_polyhedron::Vertex_const_handle& s = e->source();
    const Nef_polyhedron::Vertex_const_handle& t = e->twin()->source();
    const Nef_polyhedron::Point_3& a = s->point();
    const Nef_polyhedron::Point_3& b = t->point();
    ::glVertex3d(CGAL::to_double(a.x()),
		 CGAL::to_double(a.y()),
		 CGAL::to_double(a.z()));
    ::glVertex3d(CGAL::to_double(b.x()),
		 CGAL::to_double(b.y()),
		 CGAL::to_double(b.z()));
  }

  ::glEnd();
  GLenum gl_error = ::glGetError();
  if(gl_error != GL_NO_ERROR)
  {
    std::cerr << "OPENGL ERROR in gl_render_nef_edges!\n  "
	      << ::gluErrorString(gl_error) << "\n";
  }
}

void gl_render_nef_vertices(Nef_polyhedron* p)
{ 
  ::glBegin(GL_POINTS);
  
  for(Nef_polyhedron::Vertex_const_iterator
	v = p->vertices_begin(),
	end = p->vertices_end();
      v != end; ++v)
  {
    const Nef_polyhedron::Point_3& p = v->point();
    ::glVertex3d(CGAL::to_double(p.x()),
		 CGAL::to_double(p.y()),
		 CGAL::to_double(p.z()));
  }

  ::glEnd();
  GLenum gl_error = ::glGetError();
  if(gl_error != GL_NO_ERROR)
  {
    std::cerr << "OPENGL ERROR in gl_render_nef_vertices!\n  "
	      << ::gluErrorString(gl_error) << "\n";
  }
}
#endif // CGAL_POLYHEDRON_DEMO_USE_NEF
