// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : src/Qt_widget.C
// package       : Qt_widget
// author(s)     : Laurent Rineau & Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/Bbox_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

Qt_widget::Qt_widget(QWidget *parent, const char *name) :
  QWidget(parent, name),  Locked(0), _pointSize(4),
  _pointStyle(DISC)
{
  setCaption("CGAL::Qt_widget");

  // initialize ranges and scales
  xmin=0;
  xmax=width()-1;
  ymin=0;
  ymax=height()-1;
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  constranges=false;
  set_scales();

  // initialize the pixmap and the painter
  pixmap.resize(size());
  paint.begin(&pixmap);

  // set properties
  paint.setRasterOp(CopyROP);
  setBackgroundColor(Qt::white);
  paint.setPen(QPen(Qt::black,2));

  clear();
}

void Qt_widget::set_scales()
{
  if(!constranges)
    {
      double 
	tempmin = min(width(), height()),
	tempmax = max(xmax-xmin, ymax-ymin);
      
      xscal=yscal=(tempmin - 1)/(tempmax);
      set_scale_center(xcentre, ycentre);
    }
  else
    {
      xscal=width()/(xmax-xmin);
      yscal=height()/(ymax-ymin);
    }
}

void Qt_widget::set_scale_center(double xc, double yc)
{
  xmin = xc - (width()/xscal)/2;
  xmax = xc + (width()/xscal)/2;
  ymin = yc - (height()/yscal)/2;
  ymax = yc + (height()/yscal)/2;
  redraw();
}

void Qt_widget::resizeEvent(QResizeEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();

  paint.end();  // end painting on pixmap

  pixmap.resize(size());
  paint.begin(&pixmap); // begin again painting on pixmap

  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);

  clear();
  
  if (constranges)
    set_scales();
  else
    set_scale_center(xcentre, ycentre);
  //  emit(resized());
  redraw();
}

void Qt_widget::paintEvent(QPaintEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();


  paint.end();  // end painting on pixmap
  bitBlt(this, 0, 0, &pixmap); // copy pixmap to the Qt_widget
  paint.begin(&pixmap); // begin again painting on pixmap


  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);
}

void Qt_widget::mousePressEvent(QMouseEvent *e)
{
  emit(mousePressed(e));
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mousePressEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mousePressEvent(e);
  }
}

void Qt_widget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(mouseReleased(e));
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseReleaseEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseReleaseEvent(e);
  }
}

void Qt_widget::mouseMoveEvent(QMouseEvent *e)
{
  emit(mouseMoved(e));
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseMoveEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseMoveEvent(e);
  }
}

void Qt_widget::wheelEvent(QMouseEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->wheelEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->wheelEvent(e);
  }
}

void Qt_widget::mouseDoubleClickEvent(QMouseEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseDoubleClickEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseDoubleClickEvent(e);
  }
}

void Qt_widget::keyPressEvent(QKeyEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyPressEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyPressEvent(e);
  }
}

void Qt_widget::keyReleaseEvent(QKeyEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyReleaseEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyReleaseEvent(e);
  }
}

void Qt_widget::enterEvent(QEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->enterEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->enterEvent(e);
  }
}

void Qt_widget::leaveEvent(QEvent *e)
{
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->leaveEvent(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->leaveEvent(e);
  }
}

void Qt_widget::set_window(double  x_min, double x_max,
			   double y_min, double y_max,
			   bool const_ranges)
{
  xmin = x_min;
  xmax = x_max;
  ymin = y_min;
  ymax = y_max;
  constranges = const_ranges;
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  set_scales();
}

void Qt_widget::zoom_in(double ratio, double xc, double yc)
{
  xscal = xscal*ratio; yscal = yscal*ratio;
  xcentre = xc;
  ycentre = yc;
  set_scale_center(xcentre, ycentre);
}

void Qt_widget::zoom_in(double ratio)
{
  zoom_in(ratio,xcentre,ycentre);
}

void Qt_widget::zoom_out(double ratio)
{
  zoom_out(ratio,xcentre,ycentre);
}

void Qt_widget::zoom_out(double ratio, double xc, double yc)
{
  if(ratio!=0)
    zoom_in(1/ratio,xc,yc);
}

double Qt_widget::x_real(int x) const
{
  return(xmin+x/xscal);
}

double Qt_widget::y_real(int y) const
{
  return(ymax-y/yscal);
}

double Qt_widget::x_real_dist(double d) const
{
  return(d/xscal);
}

double Qt_widget::y_real_dist(double d) const
{
  return(d/yscal);
}

int Qt_widget::x_pixel(double x) const
{
  return( static_cast<int>((x-xmin)*xscal+0.5) );
}

int Qt_widget::y_pixel(double y) const
{
  return( - static_cast<int>((y-ymax)*yscal+0.5) );
}

int Qt_widget::x_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*xscal+0.5) );
  else
    return( static_cast<int>(d*xscal-0.5) );
}

int Qt_widget::y_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*yscal+0.5) );
  else
    return( static_cast<int>(d*yscal-0.5) );
}

Qt_widget& Qt_widget::operator<<(const Color& c)
{
  setColor(CGAL2Qt_Color(c));
  return *this;
}

Qt_widget& Qt_widget::operator<<(const PointStyle& ps)
{
  setPointStyle(ps);
  return *this;
}
	
void Qt_widget::clear() {
  painter().eraseRect(rect());
}

  Qt_widget& operator<<(Qt_widget& w, const Bbox_2& r)
  {
    int
      xmin = w.x_pixel(r.xmin()),
      ymin = w.y_pixel(r.ymin()),
      xmax = w.x_pixel(r.xmax()),
      ymax = w.y_pixel(r.ymax());

    w.painter().drawWinFocusRect(xmin, ymin, xmax-xmin, ymax-ymin);
    w.do_paint();
    return w;
  }

  void Qt_widget::attach(Qt_widget_layer *layer) {
    qt_layers.push_back(layer);
    layer->attach(this);
    layer->activate();
  }

  void Qt_widget::attach_standard(Qt_widget_layer *layer) {
    qt_standard_layers.push_back(layer);
    layer->attach(this);
    layer->activate();
  }

  bool Qt_widget::is_standard_active() {
    std::list<Qt_widget_layer*>::iterator it;
		for(it = qt_standard_layers.begin();
        it!= qt_standard_layers.end(); it++)
		  if((*it)->is_active())
        return true;
    return false;
  }


  void Qt_widget::redraw()
  {
    if(isVisible())
    {
      clear();
      lock();
        std::list<Qt_widget_layer*>::iterator it;
		    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
		      if((*it)->is_active())
			      (*it)->draw();
        for(it = qt_standard_layers.begin();
            it!= qt_standard_layers.end(); it++)
		      if((*it)->is_active())
			      (*it)->draw();
      unlock();
    }
    emit(custom_redraw());
  };
  
  void Qt_widget::detach(Qt_widget_layer* s)
  {
    qt_layers.erase(std::find(qt_layers.begin(),qt_layers.end(),s));
  }

} // namespace CGAL

// moc_source_file: ../../include/CGAL/IO/Qt_widget.h
#include "Qt_widget.moc"

#endif
