#include "MainWindow.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  MainWindow w;
  
  w.ui->viewer->restoreStateFromFile();
  
  w.show();

  a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
  a.connect(w.ui->actionExit, SIGNAL(triggered()), &a, SLOT(quit()));
  
  return a.exec();
}
