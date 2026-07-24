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
               "CGAL Mesh_3 Demo",
               "CGAL Lab (Mesh_3)",
               QStringList() << "Mesh_3");
  return app.try_exec();
}
