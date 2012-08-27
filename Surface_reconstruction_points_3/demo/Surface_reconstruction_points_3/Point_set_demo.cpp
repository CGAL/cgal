#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setOrganizationDomain("cgal.org");
  app.setOrganizationName("CGAL");
  app.setApplicationName("Point Set Demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_QT4_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  if(!args.empty() && args[0] == "--use-meta")
  {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
    args.removeAt(0);
  }

  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }
  return app.exec();
}

#ifndef USE_FORWARD_DECL
#  include "Scene.cpp"
#  include "Scene_item.cpp"
#  include "Scene_moc.cpp"
#  include "Viewer.cpp"
#  include "Viewer_moc.cpp"
#  include "MainWindow.cpp"
#  include "MainWindow_moc.cpp"
#endif
