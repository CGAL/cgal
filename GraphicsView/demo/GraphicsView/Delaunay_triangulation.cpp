

#include <QApplication>
#include "Delaunay_triangulation_MainWindow.h"


int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  Delaunay_triangulation_MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
