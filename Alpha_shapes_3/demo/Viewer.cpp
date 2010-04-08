#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>




void
Viewer::sceneChanged()
{

  Iso_cuboid_3 bb = CGAL::bounding_box(scene->points.begin(), scene->points.end());
   
  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(), bb.ymin(), bb.zmin()),
				      qglviewer::Vec(bb.xmax(),
						     bb.ymax(),
						     bb.zmax()));
    
  this->showEntireScene();
}


void
Viewer::draw()
{

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
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface();

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
Viewer::gl_draw_surface()
{
  ::glColor3f(1.0f, 0.0f, 0.0f);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_POINT_SMOOTH);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for(std::list<Point_3>::iterator it = scene->points.begin();
      it != scene->points.end();
      ++it){
    ::glVertex3d(it->x(), it->y(), it->z());
  }
  ::glEnd();
  ::glDisable(GL_POINT_SMOOTH);

  ::glEnable(GL_LIGHTING);
  ::glBegin(GL_TRIANGLES);

  ::glColor3f(0.2f, 1.0f, 0.2f);

  std::list<Facet> facets;
  scene->alpha_shape.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
  
  for(std::list<Facet>::iterator fit = facets.begin();
      fit != facets.end();
      ++fit) {
    const Cell_handle& ch = fit->first;
    const int index = fit->second;
    
    //const Vector_3& n = ch->normal(index); // must be unit vector
    
    const Point_3& a = ch->vertex((index+1)&3)->point();
    const Point_3& b = ch->vertex((index+2)&3)->point();
    const Point_3& c = ch->vertex((index+3)&3)->point();
   
    Vector_3 v = CGAL::unit_normal(a,b,c);


    ::glNormal3d(v.x(),v.y(),v.z());
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
  }

  
  ::glEnd();

}

#include "Viewer.moc"
