#include "Viewer.h"
#include <qapplication.h>

int main(int argc, char** argv)
{
  // Read command lines arguments.
  QApplication application(argc,argv);

  // Instantiate the viewer.
  Viewer viewer;
  //for Windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  application.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
  viewer.setWindowTitle("Intersection points of randomly generated circles.");

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
}
