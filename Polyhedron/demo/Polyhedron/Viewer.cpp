#include "Viewer.h"
#include <CGAL/gl.h>
#include "Scene_draw_interface.h"
#include <QMouseEvent>
#include <QKeyEvent>
#include <QGLViewer/manipulatedCameraFrame.h>
#include <QDebug>
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
  setKeyDescription(Qt::Key_T,
                    tr("Turn the camera by 180 degrees"));
  setKeyDescription(Qt::Key_M,
                    tr("Toggle macro mode: useful to view details very near from the camera, "
                       "but decrease the z-buffer precision"));
#if QGLVIEWER_VERSION >= 0x020501
  //modify mouse bindings that have been updated
  setMouseBinding(Qt::Key(0), Qt::NoModifier, Qt::LeftButton, RAP_FROM_PIXEL, true, Qt::RightButton);
  setMouseBindingDescription(Qt::ShiftModifier, Qt::RightButton,
                             tr("Select and pop context menu"));
  setMouseBinding(Qt::Key_R, Qt::NoModifier, Qt::LeftButton, RAP_FROM_PIXEL);
  //use the new API for these
  setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, SELECT);
  setMouseBindingDescription(Qt::Key(0), Qt::ShiftModifier, Qt::LeftButton,
                             tr("Selects and display context "
                                "menu of the selected item"));
#else
  setMouseBinding(Qt::SHIFT + Qt::LeftButton, SELECT);
  setMouseBindingDescription(Qt::SHIFT + Qt::RightButton,
                             tr("Selects and display context "
                                "menu of the selected item"));
#endif // QGLVIEWER_VERSION >= 2.5.0
  for(int i=0; i<16; i++)
      pickMatrix_[i]=0;
  pickMatrix_[0]=1;
  pickMatrix_[5]=1;
  pickMatrix_[10]=1;
  pickMatrix_[15]=1;
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
  glEnable(GL_DEPTH_TEST);
  d->inFastDrawing = false;
  QGLViewer::draw();
  d->draw_aux(false, this);
}

void Viewer::fastDraw()
{
  d->inFastDrawing = true;
  QGLViewer::fastDraw();
  d->draw_aux(false, this);
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  initializeOpenGLFunctions();
  glDrawArraysInstanced = (PFNGLDRAWARRAYSINSTANCEDARBPROC)this->context()->getProcAddress("glDrawArraysInstancedARB");
  if(!glDrawArraysInstanced)
  {
      qDebug()<<"glDrawArraysInstancedARB : extension not found. Spheres will be displayed as points.";
      extension_is_found = false;
  }
  else
      extension_is_found = true;

  glVertexAttribDivisor = (PFNGLVERTEXATTRIBDIVISORARBPROC)this->context()->getProcAddress("glVertexAttribDivisorARB");
  if(!glDrawArraysInstanced)
  {
      qDebug()<<"glVertexAttribDivisorARB : extension not found. Spheres will be displayed as points.";
      extension_is_found = false;
  }
  else
      extension_is_found = true;


  setBackgroundColor(::Qt::white);
  d->scene->initializeGL();

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

  viewer->glLineWidth(1.0f);
  viewer->glPointSize(2.f);
  viewer->glEnable(GL_POLYGON_OFFSET_FILL);
  viewer->glPolygonOffset(1.0f,1.0f);
  viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  viewer->glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

  if(twosides)
    viewer->glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    viewer->glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  if(antialiasing)
  {
    viewer->glEnable(GL_BLEND);
    viewer->glEnable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    viewer->glDisable(GL_BLEND);
    viewer->glDisable(GL_LINE_SMOOTH);
    viewer->glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    viewer->glBlendFunc(GL_ONE, GL_ZERO);
  }
  if(with_names)
    scene->drawWithNames(viewer);
  else
    scene->draw(viewer);
  viewer->glDisable(GL_POLYGON_OFFSET_FILL);
  viewer->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
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
    Q_EMIT selectedPoint(point.x,
                       point.y,
                       point.z);
    Q_EMIT selected(this->selectedName());
    const qglviewer::Vec orig = camera()->position();
    const qglviewer::Vec dir = point - orig;
    Q_EMIT selectionRay(orig.x, orig.y, orig.z,
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

/**
 * @brief Viewer::pickMatrix
 * Source code of gluPickMatrix slightly modified : instead of multiplying the current matrix by this value,
 * sets the viewer's pickMatrix_ so that the drawing area is only around the cursor. This is because since CGAL 4.7,
 * the drawing sustem changed to use shaders, and these need this value. pickMatrix_ is passed to the shaders in
 * Scene_item::attrib_buffers(Viewer_interface* viewer, int program_name).
 * @param x
 * @param y
 * @param width
 * @param height
 * @param viewport
 */

void Viewer::pickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
GLint viewport[4])
{
 //GLfloat m[16];
 GLfloat sx, sy;
 GLfloat tx, ty;

 sx = viewport[2] / width;
 sy = viewport[3] / height;
 tx = (viewport[2] + 2.0 * (viewport[0] - x)) / width;
 ty = (viewport[3] + 2.0 * (viewport[1] - y)) / height;

 #define M(row, col) pickMatrix_[col*4+row]
  M(0, 0) = sx;
  M(0, 1) = 0.0;
  M(0, 2) = 0.0;
  M(0, 3) = tx;
  M(1, 0) = 0.0;
  M(1, 1) = sy;
  M(1, 2) = 0.0;
  M(1, 3) = ty;
  M(2, 0) = 0.0;
  M(2, 1) = 0.0;
  M(2, 2) = 1.0;
  M(2, 3) = 0.0;
  M(3, 0) = 0.0;
  M(3, 1) = 0.0;
  M(3, 2) = 0.0;
  M(3, 3) = 1.0;
 #undef M

 //pickMatrix_[i] = m[i];
}
void Viewer::beginSelection(const QPoint &point)
{
    QGLViewer::beginSelection(point);
    //set the picking matrix to allow the picking
    static GLint viewport[4];
    camera()->getViewport(viewport);
    pickMatrix(point.x(), point.y(), selectRegionWidth(), selectRegionHeight(), viewport);

}
void Viewer::endSelection(const QPoint& point)
{
  QGLViewer::endSelection(point);
   //set dthe pick matrix to Identity
    for(int i=0; i<16; i++)
        pickMatrix_[i]=0;
    pickMatrix_[0]=1;
    pickMatrix_[5]=1;
    pickMatrix_[10]=1;
    pickMatrix_[15]=1;
}
