#include "viewer.h"

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
