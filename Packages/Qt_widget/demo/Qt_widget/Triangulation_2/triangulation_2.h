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
// file          : triangulation_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu 
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#include <fstream>
#include <stack>
#include <set>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>



#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"
#include "Qt_widget_toolbar_layers.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>	    Rep;

typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Segment_2<Rep>		    Segment;
typedef CGAL::Line_2<Rep>		    Line;
typedef CGAL::Triangle_2<Rep>		    Triangle;
typedef CGAL::Circle_2<Rep>		    Circle;

typedef CGAL::Triangulation_2<Rep>	    Triangulation;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;



typedef Delaunay::Face_handle               Face_handle;
typedef Delaunay::Vertex_handle             Vertex_handle;
typedef Delaunay::Edge                      Edge;
typedef Triangulation::Line_face_circulator Line_face_circulator;

extern int		current_state;

class Window : public QMainWindow
{
  Q_OBJECT
public:
	Window(int w, int h);
  void	set_window(double xmin, double xmax,
		    double ymin, double ymax);

public slots:
  void	new_instance();
private slots:
  void	get_new_object(CGAL::Object obj);
  void	insert_after_show_conflicts(QMouseEvent*);
  void	about();
  void	aboutQt();
  void	new_window();
  void	timerDone();
  void	generate_triangulation();
  void	save_triangulation();
  void	load_triangulation();
private:
	  void show_conflicts(Point p);
  inline  void something_changed(){current_state++;};


  CGAL::Qt_widget	  *widget;		
  CGAL::Tools_toolbar	  *newtoolbar;
  CGAL::Layers_toolbar	  *vtoolbar;
  CGAL::Standard_toolbar  *stoolbar;
  bool			  got_point;	
	  //if a CGAL::Point is received should be true
  int			  old_state;
};//endclass
