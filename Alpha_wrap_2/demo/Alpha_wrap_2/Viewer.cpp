#include <QtGui>

#include "scene.h"
#include "Viewer.h"

Viewer::Viewer(QWidget *pParent)
  : CGAL::QGLViewer(pParent)
{
  m_scene = NULL;
  m_center_x = m_center_y = 0.5;
  m_scale = 1.0;

  setAutoFillBackground(false);
  setFocusPolicy(Qt::StrongFocus);
}

void Viewer::resizeGL(int width, int height) {
  glViewport(0, 0, width, height);
  double aspect_ratio = double(height) / double(width);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.0, 1.0, -aspect_ratio, aspect_ratio, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void Viewer::initializeGL()
{
  initializeOpenGLFunctions();

  glClearColor(1., 1., 1., 0.);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_SMOOTH);
}

void Viewer::paintGL()
{
  glClear(GL_COLOR_BUFFER_BIT);
  if(!m_scene) return;

  glPushMatrix();
  glScaled(m_scale, m_scale, m_scale);
  glTranslated(-m_center_x, -m_center_y, 0.0);
  m_scene->render();
  glPopMatrix();
}

void Viewer::wheelEvent(QWheelEvent *event) {
  if(!m_scene) return;
  m_scale += (m_scale/100) * (event->angleDelta().y() / 120);
  if(m_scale < 0.0) m_scale = 0.0;

  update();
}

void Viewer::mousePressEvent(QMouseEvent *event) {
  if(!m_scene) return;
  m_mouse_click = event->pos();

  if(event->button() == Qt::LeftButton)
  {
    is_drawing = true;
    setCursor(QCursor(Qt::PointingHandCursor));
    sample_mouse_path(m_mouse_click, true,
                      Qt::ShiftModifier == QApplication::keyboardModifiers());
  }
  else
  {
    setCursor(QCursor(Qt::ClosedHandCursor));
  }
}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
  if(!m_scene) return;
  m_mouse_move = event->pos();

  if(event->buttons() == Qt::LeftButton)
  {
    if(m_mouse_move != m_mouse_click)
      sample_mouse_path(m_mouse_move, false,
                        Qt::ShiftModifier == QApplication::keyboardModifiers());
  }
  else
  {
    move_camera(m_mouse_click, m_mouse_move);
  }

  m_mouse_click = m_mouse_move;
  update();
}

void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
  if(!m_scene) return;
  m_mouse_move = event->pos();

  if(event->button() == Qt::LeftButton)
  {
    is_drawing = false;
    if(m_mouse_move != m_mouse_click)
      sample_mouse_path(m_mouse_move, true,
                        Qt::ShiftModifier == QApplication::keyboardModifiers());
  }
  else
  {
    move_camera(m_mouse_click, m_mouse_move);
  }

  m_mouse_click = m_mouse_move;
  setCursor(QCursor(Qt::ArrowCursor));
  update();
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
  if((event->key() == Qt::Key_Shift) && is_drawing) {
    m_scene->close_input();
    update();
  }
}

void Viewer::keyReleaseEvent(QKeyEvent *event)
{
  if((event->key() == Qt::Key_Shift) && is_drawing ) {
    m_scene->open_input();
    update();
  }
}

void Viewer::sample_mouse_path(const QPoint& point, bool new_cmp, bool is_closed)
{
  double x, y;
  convert_to_world_space(point, x, y);

  m_scene->add_vertex(Point_2(x, y), new_cmp, is_closed);
  m_scene->set_mouse_pos(Point_2(x, y));
}

void Viewer::move_camera(const QPoint& p0, const QPoint& p1)
{
  m_center_x -= double(p1.x() - p0.x()) / double(width());
  m_center_y += double(p1.y() - p0.y()) / double(height());
}

void Viewer::convert_to_world_space(const QPoint& point, double &x, double &y)
{
  double aspect_ratio = double(height()) / double(width());

  x = double(point.x()) / double(width());
  x = (2.0*x - 1.0) / m_scale;
  x += m_center_x;

  y = 1.0 - double(point.y()) / double(height());
  y = (2.0*y - 1.0) * aspect_ratio / m_scale;
  y += m_center_y;
}
