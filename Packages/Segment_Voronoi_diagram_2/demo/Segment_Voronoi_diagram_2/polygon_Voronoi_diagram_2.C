#ifdef CGAL_USE_QT

#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <qmainwindow.h>

//#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_segment.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_simple_polygon.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <CGAL/Timer.h>

#include <list>
#include <vector>

#include "pvd_typedefs.h"
#include "pvd_insert.h"
#include "pvd_draw.h"

#include "qt_file_toolbar.h"
#include "qt_layers_toolbar.h"
#include "qt_layers.h"


//************************************
// global variables
//************************************
SVD_2 svd;
int num_selected;
std::vector<Site> sitelist;

#include "my_window.h"

#include "qt_file_toolbar.moc"
#include "qt_layers_toolbar.moc"
#include "my_window.moc"

int
main(int argc, char* argv[])
{
  int size = 750;

  QApplication app( argc, argv );
  My_Window W(size,size);
  app.setMainWidget( &W );
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
#endif
  W.show();
  W.set_window(0,size,0,size);
  W.setCaption("Segment Voronoi diagram 2");
  return app.exec();
}


#else

#include <iostream>

int
main(int argc, char* argv[])
{
  std::cerr << "This demo needs CGAL's Qt_widget installed "
            << "in order to run..."
            << std::endl << std::endl;
  return 0;
}


#endif
