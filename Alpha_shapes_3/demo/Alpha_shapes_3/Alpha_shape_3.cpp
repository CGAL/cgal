
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>


#include <CGAL/Qt/resources.h>
#include <CGAL/Qt/init_ogl_context.h>

int main(int argc, char** argv)
{
  CGAL::Qt::init_ogl_context(4,3);

  QApplication application(argc,argv);
  application.setOrganizationDomain("geometryfactory.com");
  application.setOrganizationName("GeometryFactory");
  application.setApplicationName("Alpha Shape Reconstruction");

  // Import resources from libCGALQt (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE

  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Alpha_shape_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
