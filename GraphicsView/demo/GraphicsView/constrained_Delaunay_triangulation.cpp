

#include <QApplication>
#include "Constrained_Delaunay_triangulation_MainWindow.h"


int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  Constrained_Delaunay_triangulation_MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
