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
// file          : include/CGAL/IO/Qt_Window_tool.h
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_TOOL_H
#define CGAL_QT_WINDOW_TOOL_H

#include <CGAL/IO/Qt_Window.h>
#include <CGAL/Object.h>
#include <qobject.h>

namespace CGAL {

class Qt_widget_tool : public QObject
{
  Q_OBJECT
public:

  Qt_widget_tool();

  // attach a Qt_widget to the tool
  void attach(Qt_widget *w);
  // detach it
  void detach();

  inline bool is_attached() const
  { return (widget==0); };

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

signals:
  void new_object(Object);

protected:
  virtual void attaching()=0;
  virtual void detaching()=0;

private:
  Qt_widget *widget;
};

} // namespace CGAL

#endif // CGAL_QT_WINDOW_TOOL_H
