// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent Rineau and Radu Ursu

#include <CGAL/basic.h>


#include <CGAL/Bbox_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

Qt_widget::Qt_widget(QWidget *parent, const char *name) :
  QWidget(parent, name), set_scales_to_be_done(false), Locked(0),
  _pointSize(4), _pointStyle(DISC) 
{ 
  setCaption("CGAL::Qt_widget");

  // initialize ranges and scales
  xmin_old = xmin = -1;
  xmax_old = xmax = 1;
  ymin_old = ymin = -1;
  ymax_old = ymax = 1;
  constranges=false;
  set_scales();
  emit(rangesChanged());

  // initialize the pixmap and the painter
  painter = new QPainter;
  printer = new QPrinter;
  pixmap = new QPixmap;
  matrix = new QWMatrix;

  pixmap->resize(size());
  painter->begin(pixmap);
  painter->setWorldMatrix(*matrix);

  // set properties
  painter->setRasterOp(CopyROP);
  setBackgroundColor(Qt::white);
  painter->setPen(QPen(Qt::black,2));

  clear();
}

void Qt_widget::set_scales()
{
  if( ! isVisible() )
    {
      set_scales_to_be_done = true;
      return;
    };
  set_scales_to_be_done = false;

  if(!constranges)
    {
      xscal = yscal = min( width() / (xmax - xmin),
			   height() / (ymax - ymin) );
      double xcenter = xmin + (xmax - xmin) / 2;
      double ycenter = ymin + (ymax - ymin) / 2;

      if(xscal<1) {
        // if xscal < 1, width()/xscal > width(). then we can round it 
        // with loosing precision.
	      xmin = xcenter - (int)(width()/xscal)/2;
	      xmax = xcenter + (int)(width()/xscal)/2;
	      ymin = ycenter - (int)(height()/yscal)/2;
	      ymax = ycenter + (int)(height()/yscal)/2;
      } else {
	      xmin = xcenter - (width()/xscal)/2;
	      xmax = xcenter + (width()/xscal)/2;
	      ymin = ycenter - (height()/yscal)/2;
	      ymax = ycenter + (height()/yscal)/2;
      }
    }
    else
    {
      xscal=width()/(xmax-xmin);
      yscal=height()/(ymax-ymin);
    }
}

void Qt_widget::move_center(const double distx, const double disty)
{
  xmin += distx; xmin_old += distx;
  xmax += distx; xmax_old += distx;
  ymin += disty; ymin_old += disty;
  ymax += disty; ymax_old += disty;
  redraw();
  emit(rangesChanged());
}
void Qt_widget::set_center(const double x, const double y)
{
  if (set_scales_to_be_done) return;

  if(xscal<1) {
    xmin = x - (int)(width()/xscal)/2;
    xmax = x + (int)(width()/xscal)/2;
    ymin = y - (int)(height()/yscal)/2;
    ymax = y + (int)(height()/yscal)/2;
  } else {
    xmin = x - (width()/xscal)/2;
    xmax = x + (width()/xscal)/2;
    ymin = y - (height()/yscal)/2;
    ymax = y + (height()/yscal)/2;
  }
  xmin_old = xmin;
  xmax_old = xmax;
  ymin_old = ymin;
  ymax_old = ymax;  
  redraw();
  emit(rangesChanged());
}

void Qt_widget::resize_pixmap()
{
  // save paint state
  QFont f=painter->font();
  QBrush b=painter->brush();
  QPen p=painter->pen();
  QColor bc=painter->backgroundColor();
  QWMatrix bm = painter->worldMatrix();

  painter->end();  // end painting on pixmap
  pixmap->resize(size());
  painter->begin(pixmap); // begin again painting on pixmap
  clear();

  // restore paint state
  painter->setFont(f);
  painter->setBrush(b);
  painter->setPen(p);
  painter->setBackgroundColor(bc);
  painter->setWorldMatrix(bm);
}

void Qt_widget::resizeEvent(QResizeEvent*)
{
  resize_pixmap();
  xmin = xmin_old;
  xmax = xmax_old;
  ymin = ymin_old;
  ymax = ymax_old;
  set_scales();
  redraw();
}

void Qt_widget::showEvent(QShowEvent* e)
{
  if( set_scales_to_be_done )
    set_scales();

  return QWidget::showEvent(e);
}

void Qt_widget::paintEvent(QPaintEvent*)
{
  // save paint state
  QFont f=painter->font();
  QBrush b=painter->brush();
  QPen p=painter->pen();
  QColor bc=painter->backgroundColor();
  QWMatrix bm = painter->worldMatrix();

  painter->end();  // end painting on pixmap
  bitBlt(this, 0, 0, pixmap); // copy pixmap to the Qt_widget
  painter->begin(pixmap); // begin again painting on pixmap
  painter->setWorldMatrix(bm);

  // restore paint state
  painter->setFont(f);
  painter->setBrush(b);
  painter->setPen(p);
  painter->setBackgroundColor(bc);
}

void Qt_widget::mousePressEvent(QMouseEvent *e)
{
  emit(s_mousePressEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mousePressEvent(e);
  } 
  if(is_standard_active()){
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mousePressEvent(e);
  }
}

void Qt_widget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(s_mouseReleaseEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseReleaseEvent(e);
  }
  if(is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseReleaseEvent(e);
  }
}

void Qt_widget::mouseMoveEvent(QMouseEvent *e)
{
  emit(s_mouseMoveEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseMoveEvent(e);
  }
  if(is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseMoveEvent(e);
  }
}

void Qt_widget::wheelEvent(QWheelEvent *e)
{
  emit(s_wheelEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->wheelEvent(e);
  }
  if(is_standard_active()){
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->wheelEvent(e);
  }
}

void Qt_widget::mouseDoubleClickEvent(QMouseEvent *e)
{
  emit(s_mouseDoubleClickEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseDoubleClickEvent(e);
  }
  if(is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->mouseDoubleClickEvent(e);
  }
}

void Qt_widget::keyPressEvent(QKeyEvent *e)
{
  emit(s_keyPressEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyPressEvent(e);
  }
  if(is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyPressEvent(e);
  }
}

void Qt_widget::keyReleaseEvent(QKeyEvent *e)
{
  emit(s_keyReleaseEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyReleaseEvent(e);
  }
  if(is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->keyReleaseEvent(e);
  }
}

void Qt_widget::enterEvent(QEvent *e)
{
  emit(s_enterEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->enterEvent(e);
  }
  if(is_standard_active()){
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->enterEvent(e);
  }
}

void Qt_widget::leaveEvent(QEvent *e)
{
  emit(s_leaveEvent(e));
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->leaveEvent(e);
  }
  if(is_standard_active()){
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->leaveEvent(e);
  }
}

bool Qt_widget::event(QEvent *e)
{
  emit(s_event(e));
  QWidget::event(e);
  if(!does_standard_eat_events()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->event(e);
  }
  if(is_standard_active()){
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->event(e);
  }
  return true;
}

void Qt_widget::set_window(const double x_min, const double x_max,
			   const double y_min, const double y_max,
			   bool const_ranges)
{
  xmin_old = xmin = x_min;
  xmax_old = xmax = x_max;
  ymin_old = ymin = y_min;
  ymax_old = ymax = y_max;
  constranges = const_ranges;
  set_scales();
  redraw();
  emit(rangesChanged());
}


void Qt_widget::zoom(double ratio, double xc, double yc)
{  
  xscal = xscal*ratio; yscal = yscal*ratio;
  set_center(xc, yc);
}

void Qt_widget::zoom(double ratio)
{
  zoom(ratio,
       xmin + (xmax - xmin) / 2 ,
       ymin + (ymax - ymin) / 2 );
}

double Qt_widget::x_real(int x) const
{
  if(xscal<1)
    return(xmin+(int)(x/xscal));
  else
    return (xmin+x/xscal);
}

double Qt_widget::y_real(int y) const
{
    if(yscal<1)
      return(ymax-(int)(y/yscal));
    else
      return (ymax-y/yscal);
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
  return( static_cast<int>((x-xmin)*xscal) );
}

int Qt_widget::y_pixel(double y) const
{
  return( - static_cast<int>((y-ymax)*yscal) );
}

int Qt_widget::x_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*xscal) );
  else
    return( static_cast<int>(d*xscal) );
}

int Qt_widget::y_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*yscal) );
  else
    return( static_cast<int>(d*yscal) );
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
  painter->eraseRect(rect());
}

  Qt_widget& operator<<(Qt_widget& w, const Bbox_2& r)
  {
    int
      xmin = w.x_pixel(r.xmin()),
      ymin = w.y_pixel(r.ymin()),
      xmax = w.x_pixel(r.xmax()),
      ymax = w.y_pixel(r.ymax());

    w.get_painter().drawWinFocusRect(xmin, ymin, xmax-xmin, ymax-ymin);
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
    layer->does_eat_events = true;
  }

  bool Qt_widget::is_standard_active() {
    std::list<Qt_widget_layer*>::iterator it;
      for(it = qt_standard_layers.begin();
        it!= qt_standard_layers.end(); it++)
		  if((*it)->is_active())
        return true;
    return false;
  }
  
  bool Qt_widget::does_standard_eat_events() {
    std::list<Qt_widget_layer*>::iterator it;
      for(it = qt_standard_layers.begin();
        it!= qt_standard_layers.end(); it++)
	  if((*it)->is_active())
	    if((*it)->does_eat_events == true)
              return true;
    return false;
  }

  void Qt_widget::print_to_ps(){
    if(printer->setup(this)){
      QPainter *ptemp = new QPainter();
      ptemp = painter;
      QPainter *painter_for_printer = new QPainter(printer);
      painter = painter_for_printer;
      painter->setClipping(true);
      painter->setClipRect(rect());
      lock();
        emit(redraw_on_back());
        std::list<Qt_widget_layer*>::iterator it;
		    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
		      if((*it)->is_active())
			      (*it)->draw();
        emit(custom_redraw()); //deprecated, should use the following:
	emit(redraw_on_front());
      unlock();
      delete painter;
      painter = ptemp;
    }
  }

  void Qt_widget::redraw()
  {
    if(isVisible())
    {
      clear();
      lock();
        emit(redraw_on_back());
        std::list<Qt_widget_layer*>::iterator it;
          for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
            if((*it)->is_active())
              (*it)->draw();
        for(it = qt_standard_layers.begin();
            it!= qt_standard_layers.end(); it++)
          if((*it)->is_active())
            (*it)->draw();
        emit(custom_redraw()); //deprecated, should use the following:
        emit(redraw_on_front());
      unlock();
    }
  }
  
  void Qt_widget::detach(Qt_widget_layer* s)
  {
    qt_layers.erase(std::find(qt_layers.begin(),qt_layers.end(),s));
  }

} // namespace CGAL

// moc_source_file: ../../include/CGAL/IO/Qt_widget.h
#include "Qt_widget.moc"
