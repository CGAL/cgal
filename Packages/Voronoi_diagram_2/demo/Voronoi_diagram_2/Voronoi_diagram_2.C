#ifdef CGAL_USE_QT

#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/Timer.h>

#include <list>
#include <vector>

//#include "svd_insert.h"
//#include "svd_draw.h"

#include "typedefs.h"

//************************************
// global variables
//************************************
//VVD2 vd;
int num_selected;
std::vector<Site_2> sitelist;

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
  W.setCaption(W.title());
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
