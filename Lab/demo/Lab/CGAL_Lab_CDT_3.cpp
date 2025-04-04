#include "CGALlab.h"
#include <clocale>
#include <CGAL/Qt/resources.h>
#include <QSurfaceFormat>


/*!
 * \brief defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  CGAL_Lab app(argc, argv,
               "CGAL Lab (3D Constrained Delaunay Triangulations)",
               "Polyhedron_3 demo",
               QStringList() << "IO_surface_meshes" << "CDT_3");
  return app.try_exec();
}
