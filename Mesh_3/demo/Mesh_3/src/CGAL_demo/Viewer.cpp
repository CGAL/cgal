#include <CGAL_demo/Viewer.h>
#include <CGAL_demo/Scene_draw_interface.h>
#include <CGAL/Qt/CreateOpenGLContext.h>

Viewer::Viewer(QWidget* parent, bool antialiasing)
  : QGLViewer(CGAL::Qt::createOpenGLContext(), parent),
    scene(0),
    antialiasing(antialiasing),
    twosides(false),
    mask_(false),
    ratio_(1.)
{
  setMouseTracking(true);
}

void Viewer::setScene(Scene_draw_interface* scene)
{
  this->scene = scene;
}

void Viewer::setAntiAliasing(bool b)
{
  antialiasing = b;
  updateGL();
}

void Viewer::setTwoSides(bool b)
{
  twosides = b;
  updateGL();
}

void Viewer::setMask(bool b, double r)
{
  mask_ = b;
  ratio_ = r;
  updateGL();
}

void Viewer::draw()
{
  draw_aux(false);
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  scene->initializeGL();
  initializeOpenGLFunctions();
   setBackgroundColor(::Qt::white);
}

void Viewer::draw_aux(bool with_names)
{
  QGLViewer::draw();
  glEnable(GL_DEPTH_TEST);
  if(scene == 0)
    return;

  ::glLineWidth(1.0f);
  ::glPointSize(2.f);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f,1.0f);
  ::glClearColor(1.0f,1.0f,1.0f,0.0f);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  if(twosides)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  if(antiAliasing())
  {
    ::glEnable(GL_BLEND);
    ::glEnable(GL_LINE_SMOOTH);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    ::glDisable(GL_BLEND);
    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable(GL_POLYGON_SMOOTH_HINT);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }
  if(with_names)
    scene->drawWithNames(this);
  else
    scene->draw(this);
}

void Viewer::drawWithNames()
{
  draw_aux(true);
}

void Viewer::postSelection(const QPoint& p)
{
  Q_EMIT selected(this->selectedName());
  // do a redraw
  draw();
  Q_EMIT pointSelected(p);
}

void Viewer::postDraw()
{
  QGLViewer::postDraw();
  
  if ( mask_ )
  {
    draw_mask();
  }
}

void Viewer::draw_mask()
{
  // fill grid with transparent blue
  ::glColor4f(.4f, .4f, .4f, .7f);
  
  this->startScreenCoordinatesSystem();
  
  int width = this->width();
  int height = this->height();
  double widthF = static_cast<double>(width);
  double heightF = static_cast<double>(height);
  
  ::glDisable(GL_LIGHTING);
  ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); 
  ::glEnable(GL_BLEND);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  
  // Draws the background quad
  ::glBegin(GL_QUADS);
  
  if ( widthF > (ratio_*heightF) )
  {
    int w1 = static_cast<int>( (widthF-(heightF*ratio_)) / 2 );
    int w2 = width - w1;
    
    ::glVertex2i( 0, 0);
    ::glVertex2i( 0, height);
    ::glVertex2i( w1, height);
    ::glVertex2i( w1, 0);
    
    ::glVertex2i( w2, 0);
    ::glVertex2i( w2, height);
    ::glVertex2i( width, height);
    ::glVertex2i( width, 0);
  }
  else
  {
    int h1 = static_cast<int>( (heightF-(widthF/ratio_)) / 2 );
    int h2 = height - h1;
    
    ::glVertex2i( 0, 0);
    ::glVertex2i( 0, h1);
    ::glVertex2i( width, h1);
    ::glVertex2i( width, 0);
    
    ::glVertex2i( 0, h2);
    ::glVertex2i( 0, height);
    ::glVertex2i( width, height);
    ::glVertex2i( width, h2);
  }
  ::glEnd();
  
  ::glDisable(GL_BLEND);
  this->stopScreenCoordinatesSystem();
}
