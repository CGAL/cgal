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
#include <CGAL/IO/Qt_widget_tool.h>
#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

Qt_widget::Qt_widget(QWidget *parent, const char *name) :
  QWidget(parent, name),  Locked(0), _pointSize(4),
  _pointStyle(DISC), _has_tool(false), _has_standard_tool(false),
  current_tool(0), temp_pointer(0)
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
  if (has_tool() || has_standard_tool())
    current_tool->mousePressEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->mousePressEvent(e, *this);
}

void Qt_widget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(mouseReleased(e));
  if (has_tool() || has_standard_tool())
    current_tool->mouseReleaseEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->mouseReleaseEvent(e, *this);
}

void Qt_widget::mouseMoveEvent(QMouseEvent *e)
{
  emit(mouseMoved(e));
  if (has_tool() || has_standard_tool())
    current_tool->mouseMoveEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->mouseMoveEvent(e, *this);
}

void Qt_widget::wheelEvent(QMouseEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->wheelEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
    (*it).layer->wheelEvent(e, *this);
}

void Qt_widget::mouseDoubleClickEvent(QMouseEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->mouseDoubleClickEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->mouseDoubleClickEvent(e, *this);
}

void Qt_widget::keyPressEvent(QKeyEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->keyPressEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->keyPressEvent(e, *this);
}

void Qt_widget::keyReleaseEvent(QKeyEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->keyReleaseEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->keyReleaseEvent(e, *this);
}

void Qt_widget::enterEvent(QEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->enterEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->enterEvent(e, *this);
}

void Qt_widget::leaveEvent(QEvent *e)
{
  if (has_tool() || has_standard_tool())
    current_tool->leaveEvent(e);
  std::list<togglelayer>::iterator it;
  for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
    if((*it).active)
      (*it).layer->leaveEvent(e, *this);
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

//
//Radu Ursu coding ....
//
//-----------------------------------------------------
void Qt_widget::attach_standard(Qt_widget_tool* tool) {
  if (has_standard_tool()) {
    if(has_tool()) {
      current_tool->detach();
      emit(detached_standard_tool());
      current_tool = tool;
    } else {
      current_tool->detach();
      emit(detached_standard_tool());
      current_tool = tool;
      temp_pointer = NULL;      
    }
  } else {
    if(has_tool()) {
      temp_pointer = current_tool;
      current_tool = tool;
    } else {
      current_tool = tool;
    }
    _has_standard_tool=true;
  }
  current_tool->attach(this);
}

void Qt_widget::attach(Qt_widget_tool* tool) {
  if (has_standard_tool()) {
    current_tool->detach();
    emit(detached_standard_tool());
    if (has_tool()){
      temp_pointer = NULL;
      emit(detached_tool());
    } else
      _has_tool = true;
    current_tool = tool;
    _has_standard_tool = false;
    current_tool->attach(this);
  } else {
    if (has_tool()) {
      current_tool->detach();
      emit(detached_tool());
    } else 
      _has_tool = true;
    current_tool = tool;
    current_tool->attach(this);
  }//endif
}
void Qt_widget::detach_current_standard_tool()
{
  if (has_standard_tool())
  {
    current_tool->detach();
    _has_standard_tool = false;
    if (has_tool())
      current_tool = temp_pointer;
    else
      current_tool = NULL;
  }
};
void Qt_widget::detach_current_tool()
{
  if(has_standard_tool()){
    if(has_tool()){
      temp_pointer->detach();
      temp_pointer = NULL;
    }
  } else {
    if(has_tool()){
      current_tool->detach();
      current_tool = NULL;
    }
  }
  _has_tool = false;
};


// redraw shown scenes
// ***** Should be call when:
//    - an editable scene is changed (should be call by tools)
//    - ranges are changed (resize event)
void Qt_widget::redraw()
{
  if(isVisible())
    {
      clear();
      lock();
      std::list<togglelayer>::iterator it;
      for(it = qt_toggle_layers.begin(); it!= qt_toggle_layers.end(); it++)
	if((*it).active)
	  (*it).layer->draw(*this);
      
      unlock();
      if (has_tool() || has_standard_tool())
	current_tool->widget_repainted();
      if (has_tool() && has_standard_tool())
	temp_pointer->widget_repainted();
    }
  emit(custom_redraw());
  };
  
  // add a scene in the list of displayable scenes
  void Qt_widget::add_layer(Qt_widget_layer* s)
  {
    //qt_layers.push_back(s);
    togglelayer temp;
    temp.layer = s;
    temp.active = true;
    qt_toggle_layers.push_back(temp);
  }

  void Qt_widget::activate(Qt_widget_layer* s)
  {
    togglelayer temp;
    temp.layer = s;
    std::list<togglelayer>::iterator it = qt_toggle_layers.begin();
    bool found = false;
    while(it != qt_toggle_layers.end() && !found)
    {
      if( s == (*it).layer )
      {
	(*it).active = true;
	found = true;
      } else {
	it++;
      }
    }
  }

  void Qt_widget::deactivate(Qt_widget_layer* s)
  {
    togglelayer temp;
    temp.layer = s;
    std::list<togglelayer>::iterator it = qt_toggle_layers.begin();
    bool found = false;
    while(it != qt_toggle_layers.end() && !found)
    {
      if( s == (*it).layer )
      {
	(*it).active = false;
	found = true;
      } else {
	it++;
      }
    }
  }

  // remove a scene from the list of displayable scenes
  void Qt_widget::detach(Qt_widget_layer* s)
  {
    //qt_layers.erase(std::find(qt_layers.begin(),qt_layers.end(),s));
    togglelayer temp;
    temp.layer = s;
    std::list<togglelayer>::iterator it = qt_toggle_layers.begin();
    bool found = false;
    while(it != qt_toggle_layers.end() && !found)
    {
      if( s == (*it).layer )
      {
	qt_toggle_layers.erase(it);
	found = true;
      } else {
	it++;
      }
    }
    //qt_toggle_layers.erase(std::find(qt_toggle_layers.begin(),
//	 qt_toggle_layers.end(), temp));
  }

} // namespace CGAL

// moc_source_file: ../../include/CGAL/IO/Qt_widget.h
#include "Qt_widget.moc"

#endif
