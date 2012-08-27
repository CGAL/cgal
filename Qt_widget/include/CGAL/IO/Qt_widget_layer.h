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

#ifndef CGAL_QT_WIDGET_LAYER_H
#define CGAL_QT_WIDGET_LAYER_H

#include <CGAL/IO/Qt_widget.h>
#include <qobject.h>
#include <qcursor.h>
#include <list>


namespace CGAL {

class Qt_widget_layer : public QObject {
  Q_OBJECT
public:
  Qt_widget_layer(QObject* parent = 0, const char* name = 0) 
    : QObject(parent, name), does_eat_events(false), active(false){};

  // Event handlers
  virtual void mousePressEvent(QMouseEvent *) {} ;
  virtual void mouseReleaseEvent(QMouseEvent *) {};
  virtual void wheelEvent(QWheelEvent *) {};
  virtual void mouseDoubleClickEvent(QMouseEvent *) {};
  virtual void mouseMoveEvent(QMouseEvent *) {};
  virtual void keyPressEvent(QKeyEvent *) {};
  virtual void keyReleaseEvent(QKeyEvent *) {};
  virtual void enterEvent(QEvent *) {};
  virtual void leaveEvent(QEvent *) {};
  virtual bool event(QEvent *e) {QObject::event(e); return true;};

  bool    is_active(){return active;};	//return true if this layer is active
  bool    does_eat_events;
public slots:
  virtual void draw(){};
  void    stateChanged(int);
  void    toggle(bool);
  bool    activate(); //activate and return true if it was not active
  bool    deactivate();//deactivate and return true if it was active
signals:
  void    activated(Qt_widget_layer*);
  void    deactivated(Qt_widget_layer*);
private:
  void    attach(Qt_widget *w);//attach Qt_widget to the tool
  bool    active;	//true if this layers is active
  friend class Qt_widget;
protected:
  Qt_widget  *widget;//the pointer to the widget
  virtual void activating(){};
  virtual void deactivating(){};
};

} // namespace CGAL end

#endif // CGAL_QT_WIDGET_LAYER_H
