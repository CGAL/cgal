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
// author(s)     : Ursu Radu
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
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_euclidean_traits_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_2.h>



// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_widget.h>
//#include "Qt_widget_move_point.h"
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef CGAL::Alpha_shape_euclidean_traits_2<Rp> Gt;
typedef CGAL::Alpha_shape_vertex_base_2<Gt>	  Vb;
typedef CGAL::Triangulation_face_base_2<Gt>	  Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>	  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay;


namespace CGAL {

class Tools_toolbar : public QObject
{
	Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t);

  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);
  void was_repainted();

private:
  QToolBar      *maintoolbar;
  QToolButton   *but[10];
  Qt_widget     *widget;
  QButtonGroup  *button_group;
  int			nr_of_buttons;
	
  CGAL::Qt_widget_get_point<Rp>	    pointbut;
};//end class

};//end namespace

#endif
