
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>




int main(int argc, char** argv)
{
 QApplication application(argc,argv);
 
  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Terrain");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Terrain);
  Q_INIT_RESOURCE(CGAL);
  MainWindow mw;

  if(!application.arguments().value(1).isEmpty())
    mw.open(application.arguments().value(1));

  mw.show();

  return application.exec();
}
