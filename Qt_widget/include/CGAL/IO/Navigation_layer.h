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
// Author(s)     : Radu Ursu

#ifndef CGAL_NAVIGATION_LAYER_H
#define CGAL_NAVIGATION_LAYER_H

#include <CGAL/IO/Qt_widget_layer.h>

class Navigation_layer : public CGAL::Qt_widget_layer {
public:
  Navigation_layer(QObject *parent=0, const char* name=0)
    : Qt_widget_layer(parent, name){}
  void draw(){};
protected:
  void keyPressEvent(QKeyEvent *e){
    const double dx = widget->x_real(10) - widget->x_real(0);
    const double dy = widget->y_real(10) - widget->y_real(0);
    const double ypage = widget->y_max() - widget->y_min();
    if ( e->key() == Qt::Key_Left ){
      widget->move_center(dx, 0);
    } else if ( e->key() == Qt::Key_Right ){
      widget->move_center(-dx, 0);
    } else if ( e->key() == Qt::Key_Down ){
      widget->move_center(0, -dy);
    } else if ( e->key() == Qt::Key_Up ){
      widget->move_center(0, dy);
    } else if ( e->key() == Qt::Key_Prior ){ //PageUp
      widget->move_center(0, -ypage/2);
    } else if ( e->key() == Qt::Key_Next ){ //PageDown
      widget->move_center(0, ypage/2);
    } 
  }
  void activating(){
    widget->setFocusPolicy(QWidget::ClickFocus);
  }
};

#endif
