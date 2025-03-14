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
               "PMP demo",
               "CGAL Polygon Mesh Processing Demo",
               QStringList() << "PMP");
  return app.try_exec();
}
