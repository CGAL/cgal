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
// file          : include/CGAL/IO/Qt_widget_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_TOOLBAR_LAYERS_H
#define CGAL_QT_WIDGET_TOOLBAR_LAYERS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>
#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_constraineds.h"
#include "Qt_layer_show_points.h"


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef CGAL::Constrained_Delaunay_triangulation_2<Rp>  CDT;
typedef CDT::Constraint     Constraint;


namespace CGAL {

class Layers_toolbar : public QObject
{
	Q_OBJECT
public:
  Layers_toolbar(Qt_widget *w, QMainWindow *mw, CDT *t);
  ~Layers_toolbar()
  {
    delete showT;
    delete showP;
    delete showC;
    delete showMC;
  };
  QToolBar*	toolbar(){return maintoolbar;};

signals:
  void new_object(CGAL::Object);
	
private:
  QToolBar	*maintoolbar;
  QToolButton	*but[10];
  Qt_widget	*widget;
  QMainWindow	*window;
  //Delaunay	*dt;
  QButtonGroup  *button_group;

  int	  nr_of_buttons;
	
  CGAL::Qt_layer_show_triangulation < CDT>  *showT;  
  CGAL::Qt_layer_show_points < CDT >        *showP;  
  CGAL::Qt_layer_show_constraineds < CDT >  *showC;  
  CGAL::Qt_layer_mouse_coordinates          *showMC;
};//end class

};//end namespace

#endif
