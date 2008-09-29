
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>




int main(int argc, char** argv)
{
 QApplication application(argc,argv);
 
  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Alpha Shape Reconstruction");

  MainWindow mw;
  mw.show();

  return application.exec();
}
