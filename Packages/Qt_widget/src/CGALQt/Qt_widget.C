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
  is_the_first_time = true;

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
  if(!constranges)
    {
      double tempmin = min(width(), height());
	    double tempmax = max(xmax-xmin, ymax-ymin);
      
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
  QFont f=painter->font();
  QBrush b=painter->brush();
  QPen p=painter->pen();
  QColor bc=painter->backgroundColor();
  QWMatrix bm = painter->worldMatrix();

  painter->end();  // end painting on pixmap
  pixmap->resize(size());
  painter->begin(pixmap); // begin again painting on pixmap
  clear();
  painter->setWorldMatrix(bm);

  // restore paint state
  painter->setFont(f);
  painter->setBrush(b);
  painter->setPen(p);
  painter->setBackgroundColor(bc);

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
  emit(s_mouseReleaseEvent(e));
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
  emit(s_mouseMoveEvent(e));
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
  emit(s_wheelEvent(e));
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
  emit(s_mouseDoubleClickEvent(e));
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
  emit(s_keyPressEvent(e));
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
  emit(s_keyReleaseEvent(e));
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
  emit(s_enterEvent(e));
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
  emit(s_leaveEvent(e));
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

bool Qt_widget::event(QEvent *e)
{
  emit(s_event(e));
  QWidget::event(e);
  if(!is_standard_active()) {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
      if((*it)->is_active())
        (*it)->event(e);
  } else {
    std::list<Qt_widget_layer*>::iterator it;
    for(it = qt_standard_layers.begin();
    it!= qt_standard_layers.end(); it++)
      if((*it)->is_active())
        (*it)->event(e);
  }
  return true;
}

void Qt_widget::set_window(double  x_min, double x_max,
			   double y_min, double y_max,
			   bool const_ranges)
{
  if(!is_the_first_time)
    add_to_history();
  xmin = x_min;
  xmax = x_max;
  ymin = y_min;
  ymax = y_max;
  constranges = const_ranges;
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  set_scales();
}

void Qt_widget::set_window_p(double  x_min, double x_max,
			   double y_min, double y_max)
{
  xmin = x_min;
  xmax = x_max;
  ymin = y_min;
  ymax = y_max;  
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  set_scales();
}

void Qt_widget::zoom_in(double ratio, double xc, double yc)
{
  add_to_history(); //add the current viewport to history
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
  }

  bool Qt_widget::is_standard_active() {
    std::list<Qt_widget_layer*>::iterator it;
		for(it = qt_standard_layers.begin();
        it!= qt_standard_layers.end(); it++)
		  if((*it)->is_active())
        return true;
    return false;
  }

  void Qt_widget::print_to_ps(){
    if(printer->setup(this)){
      QPainter *ptemp = new QPainter();
      ptemp = painter;
      QPainter *painter_for_printer = new QPainter(printer);
      painter = painter_for_printer;
      lock();
        std::list<Qt_widget_layer*>::iterator it;
		    for(it = qt_layers.begin(); it!= qt_layers.end(); it++)
		      if((*it)->is_active())
			      (*it)->draw();
        emit(custom_redraw());
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
