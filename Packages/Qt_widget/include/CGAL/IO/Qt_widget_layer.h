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
// file          : include/CGAL/IO/Qt_widget_layer.h
// package       : Qt_widget
// author(s)     : Laurent Rineau & Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_LAYER_H
#define CGAL_QT_WIDGET_LAYER_H

#include <CGAL/IO/Qt_widget.h>
#include <qobject.h>
#include <list>


namespace CGAL {

class Qt_widget_layer : public QObject {
  Q_OBJECT
public:
  Qt_widget_layer() : active(false){};
  // Event handlers
  virtual void mousePressEvent(QMouseEvent *) {} ;
  virtual void mouseReleaseEvent(QMouseEvent *) {};
  virtual void wheelEvent(QMouseEvent *) {};
  virtual void mouseDoubleClickEvent(QMouseEvent *) {};
  virtual void mouseMoveEvent(QMouseEvent *) {};
  virtual void keyPressEvent(QKeyEvent *) {};
  virtual void keyReleaseEvent(QKeyEvent *) {};
  virtual void enterEvent(QEvent *) {};
  virtual void leaveEvent(QEvent *) {};
  virtual bool event(QEvent *e) {QObject::event(e); return true;};

  bool    is_active(){return active;};	//return true if this layer is active

public slots:
	virtual void draw(){};
  void    stateChanged(int);
  bool    activate(); //activate and return true if it was not active
  bool    deactivate();//deactivate and return true if it was active
signals:
  void		activated(Qt_widget_layer*);
  void		deactivated(Qt_widget_layer*);
private:
  void		attach(Qt_widget *w);//attach Qt_widget to the tool
  bool		active;	//true if this layers is active
  friend class		Qt_widget;
protected:
	Qt_widget	*widget;//the pointer to the widget
  virtual void activating(){};
  virtual void deactivating(){};
};

} // namespace CGAL end

#endif // CGAL_QT_WIDGET_LAYER_H
