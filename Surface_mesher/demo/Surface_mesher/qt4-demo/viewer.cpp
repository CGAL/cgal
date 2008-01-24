#include "viewer.h"
#include "surface.h"

void Viewer::init()
{
  setBackgroundColor(Qt::white);
  glLineStipple(5, 0xaaaa);
  glDisable(GL_LINE_STIPPLE);
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
  qglviewer::Camera new_camera = *(camera ());
  new_camera.fitBoundingBox(qglviewer::Vec(xmin, ymin, zmin),
                            qglviewer::Vec(xmax, ymax, zmax));
  camera()->interpolateTo(*new_camera.frame(), 1.);
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
