#include "MainWindow.h"
#include <QApplication>
#include <QMessageBox>
#include <CGAL/Qt/resources.h>
#include <stdexcept>

class Polyhedron_demo : public QApplication
{
public:
  Polyhedron_demo(int& argc, char **argv) : QApplication(argc, argv) {}

  bool notify(QObject* receiver, QEvent* event)
  {
    try {
      return QApplication::notify(receiver, event);
    } catch (std::exception &e) {
      // find the mainwindow to spawn an error message
      Q_FOREACH (QWidget *widget, QApplication::topLevelWidgets()) {
        if(MainWindow* mw = qobject_cast<MainWindow*>(widget)) {
          QMessageBox::critical(
            mw,
            tr("Unhandled exception"),
            e.what());
          break;
        }
      }
      QApplication::restoreOverrideCursor();
    } catch (...) {
      qFatal("Unknown exception encountered. Aborting.");
    }
    return false;
  }
};

int main(int argc, char **argv)
{
  Polyhedron_demo app(argc, argv);
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
  mainWindow.load_script(QFileInfo("autostart.js"));
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
