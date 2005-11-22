

#ifndef CGAL_QT_WIDGET_DELETEPOLYGON_H
#define CGAL_QT_WIDGET_DELETEPOLYGON_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qobject.h>
#include <qpopupmenu.h>
#include <qmessagebox.h> 
#include <qcursor.h>



#include "typedefs.h"

extern bool                                      red_active; 
//extern std::list<Polygon>                        red_pgns;
//extern std::list<Polygon>                        blue_pgns;

class Qt_widget_movepolygon_helper : public CGAL::Qt_widget_layer
{
Q_OBJECT
public:


public slots:
 
  void stateChanged(int);
};


template<class R>
class Qt_widget_deletepolygon : public Qt_widget_movepolygon_helper
{
    

public:
    
    typedef typename R::RT        RT;
    
  QCursor                    cursor;
  
 

  
  //constructor
  Qt_widget_deletepolygon(const QCursor c=QCursor(Qt::crossCursor)) :
       cursor(c) {};


private:
  QCursor oldcursor;

 

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton)
    {
        //RT x, y;
        //widget->x_real(e->x(), x);
        //widget->y_real(e->y(), y);
        //Point pt(x, y);
        //Polygon containing_p;   // the polygon that cantains point pt
       
        //Kernel _traits;
        //std::list<Polygon>*   l_of_pgns;
        //if(red_active)
        //  l_of_pgns = &red_pgns;
        //else
        //  l_of_pgns = &blue_pgns;

        // if(l_of_pgns->empty())
        // {
        //   QMessageBox::warning( widget, "There are no polygons !",
        //                                 "Generate some polygons first ");
        //  return;
        // }

        //std::list<Polygon>::iterator itr ;
        //for( itr=l_of_pgns->begin(); itr!= l_of_pgns->end(); ++itr)
        //{
        //    if (CGAL::bounded_side_2( (*itr).vertices_begin() , 
        //            (*itr).vertices_end(), pt , _traits ) ==  CGAL::ON_BOUNDED_SIDE )
        //            break;
        //
    
        //} // itr hols polygon containing pt or l_of_pgns->end()

        //if(itr==l_of_pgns->end())
        //    return; // nothing to do


        //l_of_pgns->erase(itr);
        //widget->redraw();
    }
  };

 
  
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
   
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };

};

#endif // CGAL_QT_WIDGET_DELETEPOLYGON_H
