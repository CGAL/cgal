// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Luc Habert, adapted from Qt_widget_get_circle by Radu Ursu

#ifndef CGAL_IO_QT_WIDGET_GET_CONIC_H
#define CGAL_IO_QT_WIDGET_GET_CONIC_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include<CGAL/Conic_2.h>
#include<CGAL/IO/Qt_widget_Conic_2.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif


CGAL_BEGIN_NAMESPACE

template <class R>
class Qt_widget_get_conic : public Qt_widget_layer
{
public:
  typedef typename R::Point_2   Point;
  typedef typename R::RT        RT;

  Qt_widget_get_conic(const QCursor c=QCursor(Qt::crossCursor),
		       QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), cursor(c), center_valid(false),
      axis_valid(false),conic_valid(false), istate(INPUT_CENTER) {};
  void draw(){
  }

protected:
  bool is_pure(Qt::ButtonState s){
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return 0;
    else
      return 1;
  }

  void processMouseEvent(QMouseEvent *e) {
    switch (istate) {
    case INPUT_CENTER: {
      widget->x_real(e->x(),x_center);
      widget->y_real(e->y(),y_center);
      center=Point(x_center,y_center);
      center_valid=true;
    }
      break;
    case INPUT_AXIS: {
      double x,y;
      widget->x_real(e->x(),x);
      widget->y_real(e->y(),y);
      xA=x-x_center;
      yA=y-y_center;
      axis=Segment_2<R>(Point(x_center-xA,y_center-yA),Point(x,y));
      axis_valid=true;
    }
      break;
    case INPUT_POINT: {
      double xB,yB;
      widget->x_real(e->x(),xB);
      widget->y_real(e->y(),yB);
      xB-=x_center;
      yB-=y_center;
      double xA2,yA2,xB2,yB2,xAyA,xByB,yA2_xA2;
      xA2=xA*xA; yA2=yA*yA; xB2=xB*xB;yB2=yB*yB;
      xAyA=xA*yA;xByB=xB*yB;
      yA2_xA2=yA2-xA2;

      double m11,m12,m13,m21,m22,m23;
      m11=-2*xByB*xAyA-yB2*yA2_xA2;
      m12=xB2*2*xAyA+yB2*2*xAyA;
      m13=xB2*yA2_xA2-2*xAyA*xByB;
      m21=xAyA*2*xAyA+yA2*yA2_xA2;
      m22=-2*xAyA*xA2-2*xAyA*yA2;
      m23=2*xAyA*xAyA-xA2*yA2_xA2;
      double det=xA2*m11+xAyA*m12+yA2*m13;
      if (det==0) return;

      double coeffs[6]={
        m11+m21,
        m12+m22,
        m13+m23,
        0,0,0
      };
      coeffs[3]=-2*x_center*coeffs[0]-coeffs[1]*y_center;
      coeffs[4]=-2*y_center*coeffs[2]-coeffs[1]*x_center;
      coeffs[5]=-det+coeffs[0]*x_center*x_center+coeffs[2]*y_center*y_center
        +coeffs[1]*x_center*y_center;
      double M=0;
      for (int i=0;i<6;++i) {
        double x=coeffs[i]<0?-coeffs[i]:coeffs[i];
        if (x>M) M=x;
      }
      int ncoeffs[6];
      for (int i=0;i<6;++i) {
        double y=(1000000000*coeffs[i]/M);
        ncoeffs[i]=static_cast<int>(y);
      }

      conic=Conic_2<R>(ncoeffs[0],ncoeffs[2],ncoeffs[1],
                       ncoeffs[3],ncoeffs[4],ncoeffs[5]);

      conic_valid=true;
    } 
      break;
    }
  }

  void displayConic() {
    if (center_valid) {
      *widget<<center; 
    }
    if (axis_valid) {
      *widget<<axis;
    }
    if (conic_valid) {
     *widget<<conic; 
    }
  }

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON) {
      processMouseEvent(e);
      widget->lock();
      QColor old_color = widget->color();
      *widget << CGAL::GREEN;
      displayConic();
      widget->setColor(old_color);
      widget->unlock();
    }
  }

  void mouseReleaseEvent(QMouseEvent *e) {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON) {
      processMouseEvent(e);
      widget->lock();
      QColor old_color = widget->color();
      *widget << CGAL::GREEN;
      displayConic();
      widget->setColor(old_color);
      widget->unlock();
      switch (istate) {
      case INPUT_CENTER:
        istate=INPUT_AXIS;
        break;
      case INPUT_AXIS:
        istate=INPUT_POINT;
        break;
      case INPUT_POINT:
        widget->new_object(make_object(conic));
        center_valid=false;
        conic_valid=false;
        istate=INPUT_CENTER;
        break;
      }
    }
  }

  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
        RasterOp old_raster = widget->rasterOp();
        QColor old_color = widget->color();
        widget->lock();
        widget->setRasterOp(XorROP);
        *widget << CGAL::GREEN;
        displayConic();
        widget->setRasterOp(old_raster);
        widget->setColor(old_color);
        widget->unlock();
        switch (istate) {
        case INPUT_CENTER:
          return;
        case INPUT_AXIS:
          center_valid=false;
          axis_valid=false;
          istate=INPUT_CENTER;
          break;
        case INPUT_POINT:
          axis_valid=false;
          conic_valid=false;
          istate=INPUT_AXIS;
          break;

        }
    }
  }

  void leaveEvent(QEvent *)
  {
    QColor old_color = widget->color();
    RasterOp old_raster = widget->rasterOp();//save the initial raster mode
    
    widget->lock();
    widget->setRasterOp(XorROP);
    *widget << CGAL::GREEN;
    displayConic();
    widget->unlock();
    widget->setRasterOp(old_raster);
    widget->setColor(old_color);
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    QColor old_color = widget->color();
    RasterOp old_raster = widget->rasterOp();//save the initial raster mode
    widget->setRasterOp(XorROP);
    widget->lock();
    *widget << CGAL::GREEN;
    displayConic();
    processMouseEvent(e);
    displayConic();
    widget->unlock();
    widget->setRasterOp(old_raster);
    widget->setColor(old_color);
  };
  void activating()
  {
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    istate=INPUT_CENTER;
    center_valid=false;
    axis_valid=false;
    conic_valid=false;
  };
  
  void deactivating()
  {
    widget->setFocusPolicy(oldpolicy);
    widget->setCursor(oldcursor);
  };

  QCursor cursor;
  QCursor oldcursor;

  bool center_valid,axis_valid,conic_valid;
  Conic_2<R> conic;
  Point center;
  double x_center,y_center;
  Segment_2<R> axis;
  double xA,yA; //the second point (intersection with an axis of the 
            // conic), offset wrt the center

  enum input_state {
    INPUT_CENTER,
    INPUT_AXIS,
    INPUT_POINT,
  };
  input_state istate;


  QWidget::FocusPolicy	oldpolicy;
};//end class 

CGAL_END_NAMESPACE

#endif // CGAL_IO_QT_WIDGET_GET_CONIC_H
