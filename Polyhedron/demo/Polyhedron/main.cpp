#include "MainWindow.h"
#include <QApplication>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polyhedron_3 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2); // PA: sounds weird here
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }
  return app.exec();
}

#include "MainWindow_curvature_estimation.cpp"
#include "MainWindow_subdivision_methods.cpp"
#include "MainWindow_self_intersection.cpp"
#include "MainWindow_convex_hull.cpp"
#include "MainWindow_simplify.cpp"
#include "MainWindow_kernel.cpp"
#include "MainWindow_pca.cpp"
