#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>




void
Viewer::sceneChanged()
{

  Bbox_3 bb = scene->bbox;
  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(), bb.ymin(), bb.zmin()),
				      qglviewer::Vec(bb.xmax(), bb.ymax(), bb.zmax()));
    
  this->showEntireScene();
}


void
Viewer::draw()
{
  if(frame_has_been_spun) {
    const qglviewer::Vec vec = manipulatedFrame()->orientation().axis();
    scene->setNormal(Vector_3(CGAL::to_double(vec.x),
                              CGAL::to_double(vec.y),
                              CGAL::to_double(vec.z)));
    frame_has_been_spun = false;
  }

 // define material
  float	ambient[]  =   { 0.25f,
                         0.20725f,
                         0.20725f,
                         0.922f };
  float	diffuse[]  =   { 1.0f,
                         0.829f,
                         0.829f,
                         0.922f };

  float	specular[]  = {  0.296648f,
                         0.296648f,
                         0.296648f,
                         0.522f };

  float	emission[]  = {  0.3f,
                         0.3f,
                         0.3f,
                         1.0f };
  float shininess[] = {  11.264f };

  // apply material
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);

  // anti-aliasing (if the OpenGL driver permits that)
  ::glEnable(GL_LINE_SMOOTH);

  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 
  // draw surface mesh
  bool m_view_surface = true;
  bool draw_triangles_edges = true;

  if(m_view_surface)
  {
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,3.0f);

    ::glDisable(GL_LIGHTING);
    gl_draw_vertices();

  ::glEnable(GL_LIGHTING);
    gl_draw_surface();

    gl_draw_constraints();

    if(draw_triangles_edges)
    {
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      gl_draw_surface();
    }
  }

}



void 
Viewer::gl_draw_vertices()
{
  ::glColor3f(1.0f, 0.0f, 0.0f);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_POINT_SMOOTH);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for(Finite_vertices_iterator it = scene->terrain.finite_vertices_begin();
      it != scene->terrain.finite_vertices_end();
      ++it){
    ::glVertex3d(CGAL::to_double(it->point().x()),
                 CGAL::to_double(it->point().y()),
                 CGAL::to_double(it->point().z()));
  }
  ::glEnd();
  ::glDisable(GL_POINT_SMOOTH);
}

void 
Viewer::gl_draw_surface()
{
  Converter conv;

  ::glBegin(GL_TRIANGLES);

  ::glColor3f(0.2f, 1.0f, 0.2f);

  for(Finite_faces_iterator fit = scene->terrain.finite_faces_begin();
      fit != scene->terrain.finite_faces_end();
      ++fit) {
    
    const Construction_kernel::Point_3& a = conv(fit->vertex(0)->point());
    const Construction_kernel::Point_3& b = conv(fit->vertex(1)->point());
    const Construction_kernel::Point_3& c = conv(fit->vertex(2)->point());

    Construction_kernel::Vector_3 v = CGAL::unit_normal(a,b,c);

    ::glNormal3d(v.x(),v.y(),v.z());
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
  }

  
  ::glEnd();

}

void 
Viewer::gl_draw_constraints()
{

  glDisable(GL_LIGHTING);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
  glLineWidth(1.5);
  glColor3f(8.f, 0.0f , 0.f);
  Finite_edges_iterator eit = scene->terrain.finite_edges_begin(),
    eend = scene->terrain.finite_edges_end();
  glBegin(GL_LINES);
  for(;eit != eend; ++eit){
    if(scene->terrain.is_constrained(*eit)){
      Point_3 p = eit->first->vertex(Terrain::cw(eit->second))->point();
      Point_3 q = eit->first->vertex(Terrain::ccw(eit->second))->point();
      glVertex3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
      glVertex3d(CGAL::to_double(q.x()), CGAL::to_double(q.y()), CGAL::to_double(q.z()));
    }
  }
  glEnd();
  glLineWidth(1);
}


#include "Viewer.moc"
