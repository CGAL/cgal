#include "Polyhedron_demo.h"
#include <clocale>
#include <CGAL/Qt/resources.h>
#include <QSurfaceFormat>


/*!
 * \brief defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  Polyhedron_demo app(argc, argv,
                      "Mesh_3 demo",
                      "CGAL Mesh_3 Demo",
                      QStringList() << "Mesh_3");
  return app.try_exec();
}
