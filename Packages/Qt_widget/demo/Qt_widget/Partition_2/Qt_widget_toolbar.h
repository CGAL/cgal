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
// file          : include/CGAL/IO/Qt_widget_toolbar.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_TOOLBAR_H
#define CGAL_QT_WIDGET_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Partition_traits_2.h>
#include <list>

// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_simple_polygon.h>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>
#include <qbuttongroup.h>

typedef CGAL::Cartesian<CGAL::MP_Float>     Rp;
typedef CGAL::Partition_traits_2<Rp>        Traits;
typedef Traits::Point_2                     Point_2;
typedef Traits::Polygon_2                   Cgal_Polygon;


namespace CGAL {

class Tools_toolbar : public QObject
{
  Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw);

  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);

private slots:
	
private:
  QToolBar       *maintoolbar;
  QToolButton    *but[10];
  Qt_widget      *widget;
  QButtonGroup   *button_group;
  int		 nr_of_buttons;
	
  CGAL::Qt_widget_get_simple_polygon<Cgal_Polygon>
                 getsimplebut;
};//end class

};//end namespace

#endif
