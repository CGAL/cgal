#include "Polyhedron_demo.h"
#include "MainWindow.h"
#include <QMessageBox>
#include <CGAL/Qt/resources.h>
#include <stdexcept>

struct Polyhedron_demo_impl {
  bool catch_exceptions;
  QScopedPointer<MainWindow> mainWindow;

  Polyhedron_demo_impl() : catch_exceptions(true) {}
}; // end struct Polyhedron_demo_impl

Polyhedron_demo::Polyhedron_demo(int& argc, char **argv,
                                 QString application_name,
                                 QString main_window_title)
  : QApplication(argc, argv)
  , d_ptr_is_initialized(false)
  , d_ptr(new Polyhedron_demo_impl)
{
  d_ptr_is_initialized = true;
  std::cerr.precision(17);
  std::cout.precision(17);
  std::clog.precision(17);

  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  this->setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  this->setOrganizationDomain("geometryfactory.com");
  this->setOrganizationName("GeometryFactory");
  this->setApplicationName(application_name);

  d_ptr->mainWindow.reset(new MainWindow);
  MainWindow& mainWindow = *d_ptr->mainWindow;

  mainWindow.setWindowTitle(main_window_title);
  mainWindow.show();
  QStringList args = this->arguments();
  args.removeAt(0);

  if(!args.empty() && args[0] == "--use-meta")
  {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
    args.removeAt(0);
  }
  if(!args.empty() && args[0] == "--no-try-catch")
  {
    this->do_not_catch_exceptions();
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

}

Polyhedron_demo::~Polyhedron_demo() {}

void Polyhedron_demo::do_not_catch_exceptions() {
  d_ptr->catch_exceptions = false;
}

bool Polyhedron_demo::notify(QObject* receiver, QEvent* event)
{
  if(!d_ptr_is_initialized || !d_ptr->catch_exceptions)
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

int Polyhedron_demo::try_exec()
{
  // A Qt Script may have closed the main window.
  // The following loop launch app.exec() only if the main window is visible.
  if(d_ptr->mainWindow->isVisible()) {
    return this->exec();
  }
  return 0;
}
