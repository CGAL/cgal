#ifdef CGAL_USE_QT

#include "qt_layer.h"
#include "demo_tab.h"
#include <qtabwidget.h>

/*! constructor */
Qt_layer::Qt_layer( QTabWidget * bar ) :
  myBar(bar)
{}

/*! draw - activate the current page widget draw function */
void Qt_layer::draw()
{
  // We peform downcasting from QWigdet* to Qt_widget_demo_tab*
  // , as we know that only
  // Qt_widget_demo_tab objects are stored in the tab pages.
  Qt_widget_base_tab    *w_base_p = 
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  
  TraitsType t = w_base_p->traits_type;
  
  switch ( t ) {
   case SEGMENT_TRAITS:
    {
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
       (myBar->currentPage());
     w_demo_p->lock();
     w_demo_p->draw();
     w_demo_p->unlock();
     break;
    }
   case POLYLINE_TRAITS:
    {
     Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
       (myBar->currentPage());
     w_demo_p->lock();
     w_demo_p->draw();
     w_demo_p->unlock();
     break;
    }
   case CONIC_TRAITS:
    {
     Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p = 
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
       (myBar->currentPage());
     w_demo_p->lock();
     w_demo_p->draw();
     w_demo_p->unlock();
     break;
    }
  }
  
  
}



#endif // CGAL_USE_QT
