#include "MainWindow.h"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);
  Q_FOREACH(QString filename, args) {
    mainWindow.open(filename);
  }
  return app.exec();
}

#include "MainWindow.cpp"
#include "Scene.cpp"
#include "MainWindow_subdivision_methods.cpp"
