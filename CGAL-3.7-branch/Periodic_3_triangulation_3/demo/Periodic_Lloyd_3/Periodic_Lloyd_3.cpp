#include "MainWindow.h"
#include <QApplication>

int main(int argc, char** argv)
{
  QApplication application(argc,argv);
 
  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("INRIA");
  application.setApplicationName("3D Periodic Lloyd");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Periodic_Lloyd_3);
  Q_INIT_RESOURCE(CGAL);
  MainWindow mw;
  mw.show();

  return application.exec();
}
