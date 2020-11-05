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
  Polyhedron_demo app(argc, argv,
                      "Classification demo",
                      "CGAL Classification Demo",
                      QStringList() << "Classification");
  return app.try_exec();
}
