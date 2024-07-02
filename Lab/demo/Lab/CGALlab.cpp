#include "CGALlab.h"

/*!
 * \brief defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  CGAL_Lab app(argc, argv,
                      "Polyhedron_3 demo",
                      "CGAL Lab");
  return app.try_exec();
}
