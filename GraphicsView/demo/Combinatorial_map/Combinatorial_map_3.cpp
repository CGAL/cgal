
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>




int main(int argc, char** argv)
{
 QApplication application(argc,argv);
 
  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("3D Combinatorial Map");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Combinatorial_map_3);
  Q_INIT_RESOURCE(CGAL);
  MainWindow mw;
  mw.show();

  return application.exec();
}
