#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/Qt/init_ogl_context.h>

int main(int argc, char** argv)
{

  CGAL::Qt::init_ogl_context(2, 1);

  QApplication application(argc,argv);
  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("INRIA");
  application.setApplicationName("3D Periodic Lloyd");

  // Import resources from libCGAL (Qt6).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Periodic_Lloyd_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
