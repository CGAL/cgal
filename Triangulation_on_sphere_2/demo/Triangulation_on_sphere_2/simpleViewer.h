//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr
#ifndef SIMPLE_VIEWER_H
#define SIMPLE_VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/camera.h>
#include <QKeyEvent>
#include "tools.h"


#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Cartesian_converter.h>

#include "Circular_arc_3_subsampling.h"


class Viewer : public QGLViewer
{
  typedef CGAL::Exact_spherical_kernel_3 SK;
  
  typedef Kernel::Point_3 Point_3;
  typedef std::list< std::list<Point_3> > Subsampled_arcs;

  Point_3 center_;
  double radius_;
  bool draw_balls;
  bool draw_inputs;
  double min_edge_size;
  
  Subsampled_arcs subsampled_arcs;
  std::list<Point_3> inputs;
  
  GLUquadricObj *qsphere;
  
  template <class Triangulation_on_sphere>
  void build_the_boundary(const Triangulation_on_sphere&);
  
protected :
  virtual void draw();
  virtual void init();
  virtual QString helpString() const;
  virtual void keyPressEvent(QKeyEvent *e);
public:
  Viewer(QWidget * parent) : QGLViewer(parent) {}

  template <class Triangulation_on_sphere,class Iterator>
  void open(Iterator begin, Iterator end,
            Triangulation_on_sphere& T,
            Point_3 center,
            double scale)
  {
    center_ = center;
    radius_ = scale;
    draw_balls = true;
    draw_inputs = false;
    min_edge_size = scale/100.;
    subsampled_arcs.clear();
    inputs.clear();
    std::copy(begin,end,std::back_inserter(inputs));
    build_the_boundary(T);
  }

  ~Viewer(){
    gluDeleteQuadric(qsphere);
  }
};

void Viewer::keyPressEvent(QKeyEvent *e){
  if (e->key()==Qt::Key_P){
    draw_inputs=!draw_inputs;
    updateGL();
    return;
  }

  if (e->key()==Qt::Key_S){
    draw_balls=!draw_balls;
    updateGL();
    return;
  }
  
  QGLViewer::keyPressEvent(e);
}

void Viewer::draw()
{
  glPushMatrix();
  glScalef(radius_,radius_,radius_);
  glTranslatef(-center_.x(),-center_.y(),-center_.z());
  
  if (draw_inputs){
    for (std::list<Point_3>::const_iterator 
           it=inputs.begin(), end = inputs.end();
         it !=end; ++it) {
      drawing_tools<Kernel>::draw_a_point_as_sphere(qsphere,*it,Point_3(0,0,1));
    }
  }
  
  if (draw_balls)
    drawing_tools<Kernel>::drawSphere(qsphere,center_,sqrt(radius_));
  
  glDisable(GL_LIGHT0);
  for (Subsampled_arcs::iterator 
         it_arcs=subsampled_arcs.begin(), end =subsampled_arcs.end();
       it_arcs != end; ++it_arcs){
      glLineWidth(2);
      glBegin(GL_LINE_STRIP);
      glColor3f(0,1,0);
      for (std::list<Point_3>::iterator 
             it_pt=it_arcs->begin();
           it_pt!=it_arcs->end();++it_pt) 
      {
        glVertex3f(it_pt->x(),it_pt->y(),it_pt->z());
      }
      glEnd();
  }
  glEnable(GL_LIGHT0);
  glPopMatrix();
}


template <class Triangulation_on_sphere>
void Viewer::build_the_boundary(const Triangulation_on_sphere& T)
{
  for (typename Triangulation_on_sphere::All_edges_iterator 
    it=T.all_edges_begin();it!=T.all_edges_end();++it)
  {
    if ( it->first->is_ghost() && 
         it->first->neighbor(it->second)->is_ghost() )
      continue;
    
    Point_3 source=it->first->vertex( (it->second+1)%3 )->point();
    Point_3 target=it->first->vertex( (it->second+2)%3 )->point();
    subsampled_arcs.push_front(std::list<Kernel::Point_3>());
    
    Kernel::Plane_3  plane(source,target,center_);
    Kernel::Circle_3 circle(center_,radius_,plane);
    subsample_circular_arc_3<Kernel>(circle,plane,source,target,std::back_inserter(*subsampled_arcs.begin()),min_edge_size);
  }
}


void Viewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();
 
  qsphere=gluNewQuadric();
  gluQuadricOrientation(qsphere,GLU_OUTSIDE);

  //Jane
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
  glLightfv(GL_LIGHT0,GL_POSITION,lpos);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

  ::glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
}





QString Viewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

#endif
