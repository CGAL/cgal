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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Partition_traits_2.h>

// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_segment.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>
#include <qbuttongroup.h>

typedef double                              Coord_type;
typedef CGAL::Cartesian<Coord_type>         Rp;
typedef CGAL::Exact_predicates_tag          Itag;
typedef CGAL::Triangulation_vertex_base_2<Rp>  Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Rp>
                                                Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>
                                                TDS;
typedef CGAL::Exact_predicates_tag               Itag;

typedef CGAL::Constrained_Delaunay_triangulation_2<Rp, TDS, Itag>
                                                CT;
typedef CGAL::Constrained_triangulation_plus_2<CT>      CDT;
typedef CDT::Constraint     Constraint;
typedef CGAL::Partition_traits_2<Rp>        Traits;
typedef Traits::Polygon_2                   Polygon;


namespace CGAL {
class Tools_toolbar : public QObject
{
	Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw, CDT *t);

  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);
  void was_repainted();

private slots:
  void get_new_object(CGAL::Object obj) { emit(new_object(obj)); }

private:
  QToolBar		  *maintoolbar;
  QToolButton		*but[10];
  Qt_widget		  *widget;
  QButtonGroup  *button_group;
  int			      activebutton;
  bool			    is_active;
  void			    setActiveButton(int i);
  void			    addToolButton(QToolButton *b);
  int			      nr_of_buttons;
	
  CGAL::Qt_widget_get_segment<Rp>   segmentbut;
  CGAL::Qt_widget_get_point<Rp>	    pointbut;
  CGAL::Qt_widget_get_polygon<Polygon>   polygonbut;
};//end class

};//end namespace

#endif
