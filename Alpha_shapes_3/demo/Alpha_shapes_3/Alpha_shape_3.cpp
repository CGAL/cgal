
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>


#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
 QApplication application(argc,argv);
 
  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Alpha Shape Reconstruction");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_QT4_INIT_RESOURCES;
  Q_INIT_RESOURCE(Alpha_shape_3);
  
  MainWindow mw;
  mw.show();

  return application.exec();
}
