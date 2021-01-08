#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  QSurfaceFormat fmt;
  fmt.setVersion(2, 1);
  fmt.setRenderableType(QSurfaceFormat::OpenGL);
  fmt.setProfile(QSurfaceFormat::CoreProfile);
  fmt.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(fmt);
  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
  QApplication application(argc,argv);
  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("INRIA");
  application.setApplicationName("3D Periodic Lloyd");

  // Import resources from libCGAL (QT5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Periodic_Lloyd_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
