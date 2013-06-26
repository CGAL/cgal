#include "Viewer.h"
#include <CGAL/gl.h>
#include <CGAL/check_gl_error.h>
#include "Scene_draw_interface.h"
#include <QMouseEvent>
#include <QKeyEvent>

class Viewer_impl {
public:
  Scene_draw_interface* scene;
  bool antialiasing;
  bool twosides;
  bool macro_mode;
  bool inFastDrawing;

  void draw_aux(bool with_names, Viewer*);
};

Viewer::Viewer(QWidget* parent, bool antialiasing)
  : Viewer_interface(parent)
{
  d = new Viewer_impl;
  d->scene = 0;
  d->antialiasing = antialiasing;
  d->twosides = false;
  d->macro_mode = false;
  setShortcut(EXIT_VIEWER, 0);
  setMouseBinding(Qt::SHIFT + Qt::LeftButton, SELECT);
  setMouseBindingDescription(Qt::SHIFT + Qt::RightButton, 
                             tr("Selects and display context "
                                "menu of the selected item"));
  setKeyDescription(Qt::Key_T,
                    tr("Turn the camera by 180 degrees"));
  setKeyDescription(Qt::Key_M,
                    tr("Toggle macro mode: useful to view details very near from the camera, "
                       "but decrease the z-buffer precision"));
}

Viewer::~Viewer()
{
  delete d;
}

void Viewer::setScene(Scene_draw_interface* scene)
{
  d->scene = scene;
}

bool Viewer::antiAliasing() const
{
  return d->antialiasing; 
}

void Viewer::setAntiAliasing(bool b)
{
  d->antialiasing = b;
  updateGL();
}

void Viewer::setTwoSides(bool b)
{
  d->twosides = b;
  updateGL();
}

bool Viewer::inFastDrawing() const {
  return d->inFastDrawing;
}

void Viewer::draw()
{
  d->inFastDrawing = false;
  // ::glFogf(GL_FOG_END, 2*sceneRadius());
  // ::glEnable(GL_FOG);
  QGLViewer::draw();
  d->draw_aux(false, this);
  // drawLight(GL_LIGHT0);
}

void Viewer::fastDraw()
{
  d->inFastDrawing = true;
  // ::glFogf(GL_FOG_END, 2*sceneRadius());
  // ::glEnable(GL_FOG);
  QGLViewer::fastDraw();
  d->draw_aux(false, this);
  // drawLight(GL_LIGHT0);
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  setBackgroundColor(::Qt::white);
  d->scene->initializeGL();

  // ::glFogf(GL_FOG_DENSITY, 0.05f);
  // ::glHint(GL_FOG_HINT, GL_NICEST);
  // ::glFogi(GL_FOG_MODE, GL_LINEAR);
  // static const GLfloat fogColor[] = {0.5f, 0.5f, 0.5f, 1};
  // ::glFogfv(GL_FOG_COLOR, fogColor);
}

#include <QMouseEvent>

void Viewer::mousePressEvent(QMouseEvent* event)
{
  if(event->button() == Qt::RightButton &&
     event->modifiers().testFlag(Qt::ShiftModifier)) 
  {
    select(event->pos());
    requestContextMenu(event->globalPos());
    event->accept();
  }
  else {
    QGLViewer::mousePressEvent(event);
  }
}

void Viewer::keyPressEvent(QKeyEvent* e)
{
  if(!e->modifiers()) {
    if(e->key() == Qt::Key_T) {
      turnCameraBy180Degres();
      return;
    }
    else if(e->key() == Qt::Key_M) {
      d->macro_mode = ! d->macro_mode;
      if(d->macro_mode) {
        camera()->setZNearCoefficient(0.0005f);
      } else {
        camera()->setZNearCoefficient(0.005f);
      }
      this->displayMessage(tr("Macro mode: %1").
                           arg(d->macro_mode ? tr("on") : tr("off")));
      return;
    }
  }
  //forward the event to the scene (item handling of the event)
  if (! d->scene->keyPressEvent(e) )
    QGLViewer::keyPressEvent(e);
}

void Viewer::turnCameraBy180Degres() {
  qglviewer::Camera* camera = this->camera();
  using qglviewer::ManipulatedCameraFrame;

  ManipulatedCameraFrame frame_from(*camera->frame());
  camera->setViewDirection(-camera->viewDirection());
  ManipulatedCameraFrame frame_to(*camera->frame());

  camera->setOrientation(frame_from.orientation());
  camera->interpolateTo(frame_to, 0.5f);
}

void Viewer_impl::draw_aux(bool with_names, Viewer* viewer)
{
  if(scene == 0)
    return;

  ::glLineWidth(1.0f);
  ::glPointSize(2.f);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f,1.0f);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  ::glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

  if(twosides)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  if(antialiasing)
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
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    ::glBlendFunc(GL_ONE, GL_ZERO);
  }
  if(with_names)
    scene->drawWithNames(viewer);
  else
    scene->draw(viewer);
  CGAL::check_gl_error(__FILE__, __LINE__);
}

void Viewer::drawWithNames()
{
  QGLViewer::draw();
  d->draw_aux(true, this);
}

void Viewer::postSelection(const QPoint& pixel)
{
  bool found = false;
  qglviewer::Vec point = camera()->pointUnderPixel(pixel, found);
  if(found) {
    emit selectedPoint(point.x,
                       point.y,
                       point.z);
    emit selected(this->selectedName());
    const qglviewer::Vec orig = camera()->position();
    const qglviewer::Vec dir = point - orig;
    emit selectionRay(orig.x, orig.y, orig.z,
                      dir.x, dir.y, dir.z);
  }
}

bool Viewer_interface::readFrame(QString s, qglviewer::Frame& frame)
{
  QStringList list = s.split(" ", QString::SkipEmptyParts);
  if(list.size() != 7)
    return false;
  float vec[3];
  for(int i = 0; i < 3; ++i)
  {
    bool ok;
    vec[i] = list[i].toFloat(&ok);
    if(!ok) return false;
  }
  double orient[4];
  for(int i = 0; i < 4; ++i)
  {
    bool ok;
    orient[i] = list[i + 3].toDouble(&ok);
    if(!ok) return false;
  }
  frame.setPosition(qglviewer::Vec(vec[0],
                                   vec[1],
                                   vec[2]));
  frame.setOrientation(orient[0],
                       orient[1],
                       orient[2],
                       orient[3]);
  return true;
}

QString Viewer_interface::dumpFrame(const qglviewer::Frame& frame) {
  const qglviewer::Vec pos = frame.position();
  const qglviewer::Quaternion q = frame.orientation();

  return QString("%1 %2 %3 %4 %5 %6 %7")
    .arg(pos[0])
    .arg(pos[1])
    .arg(pos[2])
    .arg(q[0])
    .arg(q[1])
    .arg(q[2])
    .arg(q[3]);
}

bool Viewer::moveCameraToCoordinates(QString s, float animation_duration) {
  qglviewer::Frame new_frame;
  if(readFrame(s, new_frame)) {
    camera()->interpolateTo(new_frame, animation_duration); 
    return true;
  }
  else
    return false;
}

QString Viewer::dumpCameraCoordinates()
{
  if(camera()->frame()) {
    return dumpFrame(*camera()->frame());
  } else {
    return QString();
  }
}
