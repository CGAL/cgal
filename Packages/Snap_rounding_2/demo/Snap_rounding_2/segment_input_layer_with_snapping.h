#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_get_segment.h>



template <class R>
class Segment_input_layer : public CGAL::Qt_widget_get_segment<R>{
private:
  typedef CGAL::Qt_widget_get_segment<R>        Base;
  typedef typename Base::Point                  Point;
  typedef typename Base::RasterOp               RasterOp;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Base::widget;
  using Base::x1;
  using Base::y1;
  using Base::x2;
  using Base::y2;
  using Base::firstpoint;
  using Base::oldcursor;
  using Base::XorROP;
#endif

  std::list<Segment_2>          *seg_list;
  //true if the user selected the first vertex
  Point                         old_point;
  QCursor                       cursor;
public:
  typedef typename R::Segment_2	Segment;
  typedef typename R::Point_2   Point_2;
  typedef typename R::FT        FT;
public:
 //constructor
  Segment_input_layer(const QCursor c=QCursor(Qt::crossCursor)) :
        CGAL::Qt_widget_get_segment<R>(), cursor(c){}
  void pass_the_structure(std::list<Segment_2>* l) {
    seg_list = l;
  }
private:
  void mousePressEvent(QMouseEvent *e)
  {
    CGAL::Qt_widget_get_segment<Rep>::mousePressEvent(e);
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if(seg_list->empty()) {
        QMessageBox::warning( widget, "There are no segments in the list!",
        "Input some segments using the left mouse button!");
        return;
      }
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      Point_2 p(x, y);
      Point_2 closest_p;  
      //this point is the closest one to the mouse coordinates
      FT min_dist;
      typename std::list<Segment_2>::const_iterator it = seg_list->begin();
      min_dist = CGAL::squared_distance(p, (*it).source());
      closest_p = (*it).source();
      
      while(it!=seg_list->end())
      {
        if (min_dist > CGAL::squared_distance(p, (*it).source())) {
          min_dist = CGAL::squared_distance(p, (*it).source());
          closest_p = (*it).source();
        }
        if (min_dist > CGAL::squared_distance(p, (*it).target())) {
          min_dist = CGAL::squared_distance(p, (*it).target());
          closest_p = (*it).target();
        }
        it++;
      }
      
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(this->XorROP);
      widget->lock();
      *widget << CGAL::GREEN << CGAL::PointSize (5) 
              << CGAL::PointStyle (CGAL::DISC);
      *widget << closest_p;
      widget->unlock();
      widget->setRasterOp(old);        
      old_point = closest_p;
      if(!firstpoint){
        x1 = closest_p.x();
        y1 = closest_p.y();
        x2 = closest_p.x();
        y2 = closest_p.y();
        firstpoint = true;
      } else {
        if(x1 != closest_p.x() || y1 != closest_p.y()) {
          widget->new_object(
            CGAL::make_object(Segment(Point_2(x1, y1), 
                                      Point_2(closest_p.x(), closest_p.y()))));
          firstpoint = false;
        }
      }
    }
  }
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    firstpoint = false;
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };
};
