#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  QApplication application(argc,argv);
 
  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("INRIA");
  application.setApplicationName("3D Periodic Lloyd");
  //for windows
  application.setAttribute(Qt::AA_UseDesktopOpenGL);

  // Import resources from libCGAL (Qt4 or QT5).
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;//New for Qt5 version !
  Q_INIT_RESOURCE(Periodic_Lloyd_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
