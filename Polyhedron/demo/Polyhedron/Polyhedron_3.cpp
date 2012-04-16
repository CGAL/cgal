#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polyhedron_3 demo");

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
#ifdef QT_SCRIPT_LIB
  if(!args.empty() && args[0] == "--debug-scripts")
  {
    mainWindow.enableScriptDebugger();
    args.removeAt(0);
  }
  mainWindow.open("autostart.js", true);
#endif
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }

  // A Qt Script may have closed the main window
  // The following loop launch app.exec() only if there is a visible
  // window.
  Q_FOREACH (QWidget *widget, QApplication::topLevelWidgets()) {
    if(widget->isVisible())
      return app.exec();
  }
  return 0;
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
