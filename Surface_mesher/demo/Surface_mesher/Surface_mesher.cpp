#include "mainwindow.h"
#include <QApplication>
#include <iostream>
#include <clocale>

int main(int argc, char** argv)
{
  QApplication application(argc,argv);

  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Surface mesher Qt5 demo");

  MainWindow w;

  if(argc>1)
    w.surface_open(argv[1]);

  w.show();
  std::setlocale(LC_ALL, "C");
  return application.exec();
  std::cerr << "Exit\n";
}
