#include "Viewer.h"
#include <qapplication.h>

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
  // Read command lines arguments.
  //for Windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
  QApplication application(argc,argv);

  // Instantiate the viewer.
  Viewer viewer;
  viewer.setWindowTitle("Intersection points of randomly generated circles.");

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
}
