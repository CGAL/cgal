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
// file          : src/Qt_Window.C
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_Widget.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/IO/Qt_Widget_tool.h>
#include <CGAL/IO/Qt_Scene.h>

namespace CGAL {

Qt_widget::Qt_widget(QWidget *parent, const char *name) :
  QWidget(parent, name), initialized(false), Locked(0), _pointSize(4),
  _pointStyle(DISC), _has_tool(false), current_tool(0)
{
  setCaption("CGAL::Qt_widget");
  initialize();
  paint.begin(&pixmap);
  setBackgroundColor(Qt::white);
  paint.setPen(QPen(Qt::black,2));
  clear();
}

void Qt_widget::initialize()
{
  xmin=0;
  xmax=width()-1;
  ymin=0;
  ymax=height()-1;
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  set_scales();
  pixmap.resize(size());
  initialized=true;
}

void Qt_widget::set_scales()
{
  double 
    tempmin = min(width(), height()),
    tempmax = max(xmax-xmin, ymax-ymin);
	
  xscal=(tempmin - 1)/(tempmax);
  yscal=(tempmin - 1)/(tempmax);
  set_scale_center(xcentre, ycentre);
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

  /*
    the only difference between an initialized Qt_widget and a
    non-initialized one is here:
    if the widget has been initialized, a resizeEvent modifies the
    scalings where as it modifies z_min(), z_max() dimensions if not.
  */
  if (!isInitialized())
    initialize();
  else
  {
    pixmap.resize(size());
    //set_scales();
  }
  paint.begin(&pixmap); // begin again painting on pixmap

  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);

  clear();
  
  set_scale_center(xcentre, ycentre);
  emit(resized());
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
  if (has_tool())
    current_tool->mousePressEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->mousePressEvent(e, *this);
  }
}

void Qt_widget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(mouseReleased(e));
  if (has_tool())
    current_tool->mouseReleaseEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->mouseReleaseEvent(e, *this);
  }
}

void Qt_widget::mouseMoveEvent(QMouseEvent *e)
{
  emit(mouseMoved(e));
  if (has_tool())
    current_tool->mouseMoveEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->mouseMoveEvent(e, *this);
  }
}

void Qt_widget::wheelEvent(QMouseEvent *e)
{
  if (has_tool())
    current_tool->wheelEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->wheelEvent(e, *this);
  }
}

void Qt_widget::mouseDoubleClickEvent(QMouseEvent *e)
{
  if (has_tool())
    current_tool->mouseDoubleClickEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->mouseDoubleClickEvent(e, *this);
  }
}

void Qt_widget::keyPressEvent(QKeyEvent *e)
{
  if (has_tool())
    current_tool->keyPressEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->keyPressEvent(e, *this);
  }
}

void Qt_widget::keyReleaseEvent(QKeyEvent *e)
{
  if (has_tool())
    current_tool->keyReleaseEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->keyReleaseEvent(e, *this);
  }
}

void Qt_widget::enterEvent(QEvent *e)
{
  if (has_tool())
    current_tool->enterEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->enterEvent(e, *this);
  }
}

void Qt_widget::leaveEvent(QEvent *e)
{
  if (has_tool())
    current_tool->leaveEvent(e);
  std::list<Qt_scene*>::iterator it;
  for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++){
    (*it)->leaveEvent(e, *this);
  }
}

void Qt_widget::init(double x_min, double x_max, double y_min)
{
  double y_max=y_min+(height()-1)*(x_max-x_min)/(width()-1);
  set_window(x_min, x_max, y_min, y_max);
}

void Qt_widget::set_window(double  x_min, double x_max, double y_min, double y_max)
{
  xmin = x_min;
  xmax = x_max;
  ymin = y_min;
  ymax = y_max;
  xcentre = xmin + (xmax - xmin)/2;
  ycentre = ymin + (ymax - ymin)/2;
  set_scales();
  initialized=true;
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

Qt_widget& Qt_widget::operator<<(Qt_widget_tool* tool)
{
  if (has_tool())
    current_tool->detach();
  current_tool=tool;
  _has_tool=true;
  tool->attach(this);
  return *this;
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

/********************************************
*Ursu Radu coding ....
*
*********************************************/
Qt_widget& Qt_widget::operator>>(Qt_widget_tool &tool)
{
  if (has_tool())
    current_tool->detach();
  current_tool=&tool;
  _has_tool=true;
  current_tool->attach(this);
  return *this;
}
void Qt_widget::detach_current_tool()
{
  if (has_tool()) 
    current_tool->detach();
  _has_tool = FALSE;
};


void Qt_widget::show_scene(Qt_scene* s)
{/*
  if(scenes_to_display.find(s)!=scenes_to_display.end())
  {
    scenes_to_display[s]=true;
    redraw();
  }*/
};


// redraw shown scenes
// ***** Should be call when:
//    - an editable scene is changed (should be call by tools)
//    - ranges are changed
void Qt_widget::redraw()
{
  if(isVisible())
    {
      clear();
      lock();
      std::list<Qt_scene*>::iterator it;
      for(it = qt_scenes.begin(); it!= qt_scenes.end(); it++)
	(*it)->draw_scene(*this);
      
      emit(redrawed());
      unlock();
      if (has_tool())
	current_tool->widget_repainted();
    }
  };
  
  // add a scene in the list of displayable scenes
  void Qt_widget::add_scene(Qt_scene* s)
  {
    qt_scenes.push_back(s);
    connect(s,SIGNAL(dying(Qt_scene*)),this,SLOT(remove_scene(Qt_scene*)));
  }

  // remove a scene from the list of displayable scenes
  void Qt_widget::remove_scene(Qt_scene* s)
  {
    qt_scenes.erase(std::find(qt_scenes.begin(),qt_scenes.end(),s));
  }

} // namespace CGAL

#include "Qt_Widget.moc"

#endif
