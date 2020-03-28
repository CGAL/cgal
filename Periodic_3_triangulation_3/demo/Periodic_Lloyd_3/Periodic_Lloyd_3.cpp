#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  QApplication application(argc,argv);

  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("INRIA");
  application.setApplicationName("3D Periodic Lloyd");
  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  application.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGAL (QT5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Periodic_Lloyd_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
