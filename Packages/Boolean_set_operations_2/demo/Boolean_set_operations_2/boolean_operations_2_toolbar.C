

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>

#include "boolean_operations_2_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>
#include <CGAL/IO/pixmaps/notool.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include "demo_pointlocation.xpm"



#include <qiconset.h>


class MyWindow;
Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, 
                             QMainWindow *mw) : 
    QToolBar(mw, "NT")
  {

   
    w->attach(&getsimplebut);
    w->attach(&getcirclebut);
    w->attach(&locatebut);
    /*w->attach(delete_red_but);
    w->attach(delete_blue_but);*/
    //w->attach(&delete_polygon);
    
    getsimplebut.deactivate();
    getcirclebut.deactivate();
    locatebut.deactivate();
    /*delete_red_but->deactivate();
    delete_blue_but->deactivate();*/
    //delete_polygon.deactivate();
    //set the widget
    widget = w;

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)movepoint_small_xpm ),
                  QPixmap( (const char**)movepoint_xpm ));
    QIconSet set3(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));
    QIconSet set4(QPixmap( (const char**)pointlocation_xpm ),
                  QPixmap( (const char**)pointlocation_xpm ));
    QIconSet set5(QPixmap( (const char**)circle_xpm ),
                  QPixmap( (const char**)circle_small_xpm ));
    



  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "polygontool");
  but[1]->setIconSet(set3);
  but[1]->setTextLabel("Input Polygon");
  but[2] = new QToolButton(this, "circletool");
  but[2]->setIconSet(set5);
  but[2]->setTextLabel("Input Circle");
  but[3] = new QToolButton(this, "locatetool");
  but[3]->setIconSet(set4);
  but[3]->setTextLabel("Locate Polygon");
 
   
  int nr_of_buttons = 4;
  button_group = new QButtonGroup(0, "My_group");
  for(int i = 0; i < nr_of_buttons; i++)
  {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  
 
  connect(but[1], SIGNAL(stateChanged(int)),
        &getsimplebut, SLOT(stateChanged(int)));

  connect(but[2], SIGNAL(stateChanged(int)),
        &getcirclebut, SLOT(stateChanged(int)));

  connect(but[3], SIGNAL(stateChanged(int)),
        &locatebut, SLOT(stateChanged(int)));
}



  
#include "boolean_operations_2_toolbar.moc"

#endif
