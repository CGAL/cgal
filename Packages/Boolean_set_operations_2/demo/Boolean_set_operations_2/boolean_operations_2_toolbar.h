
#ifndef MIN_CIRCLE_2_TOOLBAR_H
#define MIN_CIRCLE_2_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>



#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_point.h>
//#include <CGAL/IO/Qt_widget_get_simple_polygon.h>
#include "Qt_widget_get_circ_polygon.h"

#include <CGAL/IO/Qt_widget_get_circle.h> 
#include <CGAL/IO/pixmaps/circle.xpm>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>
#include <qcursor.h>
#include <qradiobutton.h> 
#include <qvbuttongroup.h> 



#include "typedefs.h"
#include "Qt_widget_locate_layer.h"

extern bool                                      red_active; 
extern Polygon_set                               red_set;
extern Polygon_set                               blue_set;

class Tools_toolbar : public QToolBar
{
  Q_OBJECT
public:
	Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw);
  ~Tools_toolbar(){};

  void deactivate()
  {
    getsimplebut.deactivate();
    getcirclebut.deactivate();
    but[0]->toggle(); // toggle the 'deactivate layer button'
  }

  void reset()
  {
    locatebut.reset();
  }

private:
  QToolButton     *but[10];
  QButtonGroup    *button_group;
  CGAL::Qt_widget *widget;
  
  CGAL::Qt_widget_get_circ_polygon<Kernel>     getsimplebut;
  CGAL::Qt_widget_get_circle<Kernel>           getcirclebut;
  Qt_widget_locate_layer                       locatebut;
  
};//end class


#endif
