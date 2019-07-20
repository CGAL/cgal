#include "viewer.h"
#include "surface.h"
#include <QAction>
#include <CGAL/Qt/manipulatedCameraFrame.h>

Viewer::Viewer(QWidget* parent)
  : CGAL::QGLViewer(parent), surface(0)
{
  // Do not store state in a file
  setStateFileName("");
}

void Viewer::init()
{
  setBackgroundColor(Qt::white);
  glLineStipple(5, 0xaaaa);
  glDisable(GL_LINE_STIPPLE);

  // anti-aliasing
  glEnable(GL_BLEND);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

QString Viewer::helpString() const
{
  return ""
  "<h1>Surface mesher demo</h1>\n"
  "<p>No help availlable for now.</p>";
}

void Viewer::interpolateToFitBoundingBox(double xmin, double ymin, double zmin,
                                         double xmax, double ymax, double zmax)
{
  QAction* auto_resize = parent()->parent()->findChild<QAction*>("actionAuto_resize");
  Q_ASSERT_X(auto_resize, "Viewer::interpolateToFitBoundingBox", "cannot find action \"actionAuto_resize\"");
  if(auto_resize && auto_resize->isChecked())
  {
    CGAL::qglviewer::Camera new_camera = *(camera ());
    new_camera.fitBoundingBox(CGAL::qglviewer::Vec(xmin, ymin, zmin),
                              CGAL::qglviewer::Vec(xmax, ymax, zmax));
    camera()->interpolateTo(*new_camera.frame(), 1.);
  }
}

void Viewer::draw()
{
  if(surface)
    surface->draw();
}

void Viewer::drawWithNames()
{
  if(surface)
    surface->drawWithNames();
}
void Viewer::postSelection(const QPoint& p)
{
  if(surface)
    surface->postSelection(p);
}

void Viewer::set_surface(Surface* s)
{
  surface = s;
}

