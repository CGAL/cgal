#include "Polyhedron_demo.h"
#include <clocale>
#include <CGAL/Qt/resources.h>
#include <QSurfaceFormat>


/*!
 * \brief Defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  QSurfaceFormat fmt;

  fmt.setVersion(4, 3);
  fmt.setRenderableType(QSurfaceFormat::OpenGL);
  fmt.setProfile(QSurfaceFormat::CoreProfile);
  fmt.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(fmt);
  QStringList keywords;
  keywords << "PMP";
    Polyhedron_demo app(argc, argv,
                      "PMP demo",
                      "CGAL Polygon Mesh Processing Demo",
                        keywords);
  //We set the locale to avoid any trouble with VTK
  std::setlocale(LC_ALL, "C");
  return app.try_exec();
}
