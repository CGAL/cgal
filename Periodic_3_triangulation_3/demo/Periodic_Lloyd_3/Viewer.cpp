#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>


void
Viewer::init()
{
  setBackgroundColor(::Qt::white);
  this->camera()->setSceneBoundingBox(
      qglviewer::Vec(-1.,-1.,-1.),
      qglviewer::Vec( 1., 1., 1.));
}


void
Viewer::sceneChanged()
{
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
  ::glColor3f(1.0f, 0.72f, 0.06f);
  ::glDisable(GL_LIGHTING);

  ::glEnable(GL_BLEND);
  ::glEnable(GL_POINT_SMOOTH);
  ::glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  ::glEnable(GL_LINE_SMOOTH);
  ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for(Periodic_point_iterator ppit
	= scene->periodic_triangulation.periodic_points_begin(
	    P3DT3::UNIQUE) ;
      ppit != scene->periodic_triangulation.periodic_points_end(P3DT3::UNIQUE);
      ++ppit){
    Point_3 p(scene->periodic_triangulation.point(*ppit));
    ::glVertex3d(p.x(), p.y(), p.z());
  }
  ::glEnd();

  ::glBegin(GL_LINES);

  ::glColor3f(0.27f, 0.51f, 0.7f);

  for (Periodic_triangle_iterator ptit
	 = scene->periodic_triangulation.periodic_triangles_begin(
	     P3DT3::UNIQUE);
       ptit != scene->periodic_triangulation.periodic_triangles_end(
	   P3DT3::UNIQUE);
      ++ptit) {
    for (int i=0 ; i<4 ; i++) {
      Segment_3 dual = scene->periodic_triangulation.segment(
	  scene->periodic_triangulation.dual(*(ptit.get_facet())));

      FT sz = dual.source().z();
      FT tz = dual.target().z();

      if (scene->two_dimensional && ((sz-tz > 0.5) || (sz-tz < -0.5))) continue;
      
      if (scene->two_dimensional) { sz = 0.; tz = 0.; }
      FT sx = dual.source().x();
      FT tx = dual.target().x();
      FT sy = dual.source().y();
      FT ty = dual.target().y();

      ::glVertex3d(sx,sy,sz); ::glVertex3d(tx,ty,tz);
      Iso_cuboid_3 d = scene->periodic_triangulation.domain();
      if (scene->eight_copies) {
	::glColor3f(0.69f, 0.77f, 0.87f);
	::glVertex3d(sx+0.,sy+d.ymax()-d.ymin(),sz+0.);
	::glVertex3d(tx+0.,ty+d.ymax()-d.ymin(),tz+0.);
	::glVertex3d(sx+d.xmax()-d.xmin(),sy+0.,sz+0.);
	::glVertex3d(tx+d.xmax()-d.xmin(),ty+0.,tz+0.);
	::glVertex3d(sx+d.xmax()-d.xmin(),sy+d.ymax()-d.ymin(),sz+0.);
	::glVertex3d(tx+d.xmax()-d.xmin(),ty+d.ymax()-d.ymin(),tz+0.);
	if (!scene->two_dimensional) {
	  ::glVertex3d(sx+0.,sy+0.,sz+d.zmax()-d.zmin());
	  ::glVertex3d(tx+0.,ty+0.,tz+d.zmax()-d.zmin());
	  ::glVertex3d(sx+0.,sy+d.ymax()-d.ymin(),sz+d.zmax()-d.zmin());
	  ::glVertex3d(tx+0.,ty+d.ymax()-d.ymin(),tz+d.zmax()-d.zmin());
	  ::glVertex3d(sx+d.xmax()-d.xmin(),sy+0.,sz+d.zmax()-d.zmin());
	  ::glVertex3d(tx+d.xmax()-d.xmin(),ty+0.,tz+d.zmax()-d.zmin());
	  ::glVertex3d(sx+d.xmax()-d.xmin(),sy+d.ymax()-d.ymin(),sz+d.zmax()-d.zmin());
	  ::glVertex3d(tx+d.xmax()-d.xmin(),ty+d.ymax()-d.ymin(),tz+d.zmax()-d.zmin());
	}
	::glColor3f(0.27f, 0.51f, 0.7f);

      }
    }
  }

  ::glEnd();
  ::glEnable(GL_LIGHTING);
}

Viewer::Vec
Viewer::next_around_circle(const float& phi, const Vec& pos, const Vec& ori) {
  Vec cam = pos-ori;
  Vec cam_norm = cam/cam.norm();

  Vec y(cam_norm.z, 0, -cam_norm.x);
  Vec y_norm = y/y.norm();
  
  Vec new_cam = ori + (cam_norm*cos(phi) + y_norm*sin(phi)) * cam.norm();
  return new_cam;  
}

void
Viewer::render_video()
{
  setSnapshotFormat("PNG");
  for (int alpha=0 ; alpha <= 100 ; alpha++ ) {
    emit (valueChanged(alpha));
    std::cout<<alpha<<std::endl;
    QString alpha_str;
    alpha_str.setNum(alpha);
    displayMessage(QString("alpha: ") + alpha_str + QString("%"),10000);
    
    for (int fr=0 ; fr < 50 ; fr++) {
      Vec cam = camera()->position();
      Vec ori = sceneCenter();
      Vec new_cam = next_around_circle(0.01,cam,ori);
      camera()->setPosition(new_cam);
      camera()->lookAt(ori);
      this->showEntireScene();
      saveSnapshot(true);
    }
  }
}


#include "Viewer.moc"
