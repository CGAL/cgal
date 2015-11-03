#include "MainWindow.h"
#include <QApplication>
#include <QMessageBox>
#include <CGAL/Qt/resources.h>
#include <stdexcept>
/*!
 * \brief The Polyhedron_demo class : defines the main function
 */

class Polyhedron_demo : public QApplication
{
  bool catch_exceptions;
public:
  /*!
   * Constructor : calls the constructor of QApplication.
   */
  Polyhedron_demo(int& argc, char **argv)
    : QApplication(argc, argv)
    , catch_exceptions(true)
  {}

  void do_not_catch_exceptions() {
    catch_exceptions = false;
  }

  /*!
   * Catches unhandled exceptions from all the widgets.
   */
  bool notify(QObject* receiver, QEvent* event)
  {
    if(!catch_exceptions)
      return QApplication::notify(receiver, event);
    else try {
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

/*!
 * \brief Defines the entry point of the demo.
 * Creates the application and sets a main window.
 */
int main(int argc, char **argv)
{
  Polyhedron_demo app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Polyhedron_3 demo");
  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  app.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  if(!args.empty() && args[0] == "--use-meta")
  {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
    args.removeAt(0);
  }
  if(!args.empty() && args[0] == "--no-try-catch")
  {
    app.do_not_catch_exceptions();
    args.removeAt(0);
  }
#ifdef QT_SCRIPT_LIB
  if(!args.empty() && args[0] == "--debug-scripts")
  {
    mainWindow.enableScriptDebugger();
    args.removeAt(0);
  }
  QFileInfo autostart_js("autostart.js");
  if(autostart_js.exists()) {
    mainWindow.load_script(autostart_js);
  }
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
#  include "Viewer.cpp"
#  include "Viewer_moc.cpp"
#  include "MainWindow.cpp"
#  include "MainWindow_moc.cpp"
#endif
