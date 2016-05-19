#include "Polyhedron_demo.h"
#include <clocale>

/*!
 * \brief Defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  Polyhedron_demo app(argc, argv,
                      "Polyhedron_3 demo",
                      "CGAL Polyhedron Demo");
  //We set the locale to avoid any trouble with VTK
  std::setlocale(LC_ALL, "C");
  return app.try_exec();
}
