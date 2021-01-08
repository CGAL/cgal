#include "MainWindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
  QSurfaceFormat fmt;
  fmt.setVersion(2, 1);
  fmt.setRenderableType(QSurfaceFormat::OpenGL);
  fmt.setProfile(QSurfaceFormat::CoreProfile);
  fmt.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(fmt);
  QApplication a(argc, argv);
  MainWindow w;
  //w.ui->setupUi(w);
  w.ui->viewer->restoreStateFromFile();

  w.show();

  a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
  a.connect(w.ui->actionExit, SIGNAL(triggered()), &a, SLOT(quit()));

  return a.exec();
}
