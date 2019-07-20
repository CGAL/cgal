#include "Polyhedron_demo.h"

/*!
 * \brief Defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  Polyhedron_demo app(argc, argv,
                      "Polyhedron_3 demo",
                      "CGAL Polyhedron Demo");
  return app.try_exec();
}
