#ifndef APOLLONIUS_GRAPH_2_EDIT_VERTEX_H
#define APOLLONIUS_GRAPH_2_EDIT_VERTEX_H

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>


#include <qobject.h>
#include <qpopupmenu.h>
#include <qcursor.h>


class Edit_vertex_layer_helper
  : public CGAL::Qt_widget_layer
{
  Q_OBJECT
public:
  virtual void delete_vertexi(){};
  virtual void move_vertexi(){};
  virtual void change_weighti(){};

protected:
  void emit_signal() { 
    emit( apollonius_graph_changed() );
  }

public slots:
  void stateChanged(int i)
  {
    if( i == 2 ) {
      activate();
    } else if ( i == 0 ) {
      deactivate();
    }
  }

  void delete_vertex()
  {
    delete_vertexi();
    emit( apollonius_graph_changed() );
  }

  void move_vertex() {
    move_vertexi();
  }

  void change_weight() {
    change_weighti();
  }

signals:
  void apollonius_graph_changed();
};

#include "edit_vertex_layer.moc"


template <class AG>
class Edit_vertex_layer : public Edit_vertex_layer_helper
{
public:

  typedef typename AG::Site_2		        Site_2;
  typedef typename AG::Point_2                  Point_2;
  typedef typename AG::Face_handle              Face_handle;
  typedef typename AG::Vertex_handle            Vertex_handle;
  typedef typename AG::Geom_traits              GT;
  typedef typename GT::FT                       FT;


protected:
  FT                                                  first_x, first_y;
  FT                                                  x2, y2;
  bool                                                wasrepainted;
  bool                                                on_first; 
           //true if right mouse button was pressed
  bool
  move_button_pressed; 
           //true if the popup's move button was pressed
  bool
  change_weight_pressed; 
          //true if the popup's change_weight button was pressed
  Vertex_handle                                       current_v;
          //the vertex that will be processed
  Point_2                                             old_point;
  FT                                                  old_weight;
          //contains the old vertex that should be removed
  AG*                                                  ag;
          //pointer to regular triangulation being used
  QPopupMenu*                                          popup;
          //the popup being displayed when right mouse button is pressed
public:
  Edit_vertex_layer(AG* ag)
    : wasrepainted(true), on_first(false),
      move_button_pressed(false), change_weight_pressed(false),
      ag(ag) {};

  //  void set_apollonius_graph (AG* ag) { this->ag = ag; }
  
  template < class TRIANGULATION > 
  Vertex_handle
  closest_vertex(const TRIANGULATION &T,
		 Face_handle f,	const Point_2& p)
  {
    return ag->nearest_neighbor(p, f->vertex());
#if 0
    Vertex_handle v ;
    typename GT::Compare_distance_2 cmp =
      T.geom_traits().compare_distance_2_object();

    if( T.is_infinite(f)){
      int i = f->index(T.infinite_vertex());
      Bare_point pcwi = f->vertex(f->cw(i))->point();
      Bare_point pccwi = f->vertex(f->ccw(i))->point();
      v =  cmp(p, pcwi, pccwi) == CGAL::SMALLER ? f->vertex(f->cw(i)) :
                                                  f->vertex(f->ccw(i));
    }
    else{ 
      v = f->vertex(0);
      if (cmp(p, f->vertex(1)->point(), v->point()) == CGAL::SMALLER) 
        v = f->vertex(1);
      if (cmp(p, f->vertex(2)->point(), v->point()) == CGAL::SMALLER) 
        v = f->vertex(2);
    }
    return v;
#endif
  }

private:
  QCursor oldcursor;

  void draw() {
    wasrepainted = true;
  }

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton && on_first) {
      on_first = false;
    }

    if(e->button() == Qt::RightButton) {
      if ( ag->dimension() < 0 ) { return; }
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      Point_2 p(x, y);
      Vertex_handle v = ag->nearest_neighbor(p);
      //save the initial raster mode
      //      RasterOp old = widget->rasterOp();
      CGAL::PointStyle pstyle = widget->pointStyle();
      int psize = widget->pointSize();
      //      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::RED << CGAL::PointSize(10)
	      << CGAL::PointStyle(CGAL::CIRCLE);
      if( !wasrepainted ) {
	*widget << old_point;
	*widget << CGAL::RED;
	*widget << old_point;
      }
      *widget << v->site();
      widget->unlock();
      //      widget->setRasterOp(old);
      *widget << CGAL::PointSize(psize) << pstyle;
      popup->popup(widget->mapToGlobal(e->pos()));
      old_point = v->site().point();
      //      old_weight = v->site().weight();
      current_v = v;
      wasrepainted = false;
      on_first = false;
    }
  }
  
  void mouseMoveEvent(QMouseEvent *e)
  {
    if ( on_first ) {
      if( move_button_pressed ){
        FT x, y;
        widget->x_real(e->x(), x);
        widget->y_real(e->y(), y);
    		
        if( !wasrepainted ) {
	  *widget << old_point;
	}
        *widget << Point_2(x, y);
	FT wght = current_v->site().weight();
        ag->remove(current_v);
        current_v = ag->insert(Site_2(Point_2(x, y), wght/*old_weight*/));
        widget->redraw();	//redraw the scenes
	old_point = Point_2(x, y);
      } else if( change_weight_pressed ) {
        FT x, y;
        widget->x_real(e->x(), x);
        widget->y_real(e->y(), y);

	Point_2 lastp = current_v->site().point();
	ag->remove(current_v);
	FT w = CGAL::sqrt(CGAL::squared_distance(lastp,
						 Point_2(x,y)));
	Site_2 s(lastp, w);
        current_v = ag->insert(s);
        widget->redraw();	//redraw the scenes
	old_point = lastp;
      }
    }
  }

  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
    wasrepainted = false;
    popup = new QPopupMenu( widget, 0);
    popup->insertItem("Delete Vertex", this, SLOT(delete_vertex()));
    popup->insertItem("Move Vertex", this,  SLOT(move_vertex()));
    popup->insertItem("Change Weight", this, SLOT(change_weight()));
  }

  void deactivating()
  {
    widget->setCursor(oldcursor);
  }
  
  void delete_vertexi()
  {
    ag->remove(current_v);
    widget->redraw();	//redraw the scenes
  }
  
  void move_vertexi()
  {
    on_first = true;
    change_weight_pressed = false;
    move_button_pressed = true;
    double x = CGAL::to_double( old_point.x() );
    double y = CGAL::to_double( old_point.y() );
    widget->cursor().setPos(widget->mapToGlobal(
                            QPoint( widget->x_pixel(x), 
				    widget->y_pixel(y) )
			    ));
  }
  
  void change_weighti()
  {
    on_first = true;
    move_button_pressed = false;
    change_weight_pressed = true;
    double x = CGAL::to_double( old_point.x() );
    double y = CGAL::to_double( old_point.y() );
    widget->cursor().setPos(widget->mapToGlobal(
                            QPoint( widget->x_pixel(x), 
				    widget->y_pixel(y) )
			    ));
  }
}; //end class 




#endif // CGAL_USE_QT

#endif // APOLLONIUS_GRAPH_2_EDIT_VERTEX_H
