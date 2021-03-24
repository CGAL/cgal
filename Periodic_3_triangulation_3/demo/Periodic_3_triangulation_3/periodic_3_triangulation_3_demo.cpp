#include "MainWindow.h"

#include <QApplication>
#include <CGAL/Qt/init_ogl_context.h>

int main(int argc, char *argv[])
{
  CGAL::Qt::init_ogl_context(2,1);
  QApplication a(argc, argv);
  MainWindow w;

  w.show();

  a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
  a.connect(w.ui->actionExit, SIGNAL(triggered()), &a, SLOT(quit()));

  return a.exec();
}
