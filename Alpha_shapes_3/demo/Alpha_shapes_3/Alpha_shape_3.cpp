
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>


#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  QSurfaceFormat fmt;
#ifdef Q_OS_MAC
  fmt.setVersion(4, 1);
#else
  fmt.setVersion(4, 3);
#endif
  fmt.setRenderableType(QSurfaceFormat::OpenGL);
  fmt.setProfile(QSurfaceFormat::CoreProfile);
  fmt.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(fmt);

  //for Windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  QApplication application(argc,argv);

  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Alpha Shape Reconstruction");

  // Import resources from libCGALQt (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE

  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Alpha_shape_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
