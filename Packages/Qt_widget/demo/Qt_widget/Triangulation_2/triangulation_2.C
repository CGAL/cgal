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
// file          : triangulation_2.C
// package       : Qt_widget
// author(s)     : Radu Ursu 
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>


int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

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

const QString my_title_string("Triangulation Demo with"
			      " CGAL Qt_widget");
Delaunay	tr1;
int		current_state;

class Window : public QMainWindow
{
  Q_OBJECT
public:
  Window(int w, int h)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
    connect(widget, SIGNAL(s_mouseReleaseEvent(QMouseEvent*)), this,
          SLOT(insert_after_show_conflicts(QMouseEvent*)));    
	
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
    this, SLOT(timerDone()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Triangulation", this, 
		      SLOT(load_triangulation()), CTRL+Key_L);
    file->insertItem("&Save Triangulation", this,
		      SLOT(save_triangulation()), CTRL+Key_T);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()));
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp,
		      SLOT( closeAllWindows() ), CTRL+Key_Q );


    // drawidgetg menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate_triangulation", this,
		      SLOT(generate_triangulation()), CTRL+Key_G );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    newtoolbar = new CGAL::Tools_toolbar(widget, this, &tr1);	
    //the new scenes toolbar
    vtoolbar = new CGAL::Layers_toolbar(widget, this, &tr1);
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
    this->addToolBar(vtoolbar->toolbar(), Top, FALSE);
  
    *widget << CGAL::BackgroundColor (CGAL::BLACK);

    resize(w,h);
    widget->show();

    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    got_point = FALSE;
    old_state = 0;
  };
  void set_window(double xmin, double xmax,
			  double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }


private slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    widget->clear_history();
    tr1.clear();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }
	
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    Segment s;
    Line l;
    if (CGAL::assign(l,obj))
    {
      if (tr1.dimension()<2) return;
      widget->redraw();
      widget->lock();
      Line_face_circulator lfc = 
	tr1.line_walk(l.point(1), l.point(2)), done(lfc);
      if(lfc == (CGAL_NULL_TYPE) NULL){
      } else {
	*widget << CGAL::BLUE;
	*widget << CGAL::FillColor(CGAL::WHITE);
	do{
	  if(! tr1.is_infinite( lfc  ))
	    *widget << tr1.triangle( lfc );
	}while(++lfc != done);
      }
      *widget << CGAL::GREEN << l ;
      *widget << CGAL::noFill;
      widget->unlock();
    } else if(CGAL::assign(p,obj)) {
      got_point = TRUE;
      show_conflicts(p);
      tr1.insert(p);
    } 
  }

  void insert_after_show_conflicts(QMouseEvent*)
  {
    if(got_point)
    {
      got_point = FALSE;
      widget->redraw();
      something_changed();
    }
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Triangulation,\n"
  		"Copyright CGAL @2001");
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    Window *ed = new Window(500, 500);
    ed->setCaption("Layer");
    ed->show();
    ed->set_window(-1.1, 1.1, -1.1, 1.1);
    something_changed();
  }

  void timerDone()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void generate_triangulation()
  {
    tr1.clear();
    widget->clear_history();
    widget->lock();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); // set the Visible Area to the Interval

    // send resizeEvent only on show.
    widget->unlock();
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++)
      tr1.insert(*g++);  
    widget->redraw();
    something_changed();
  }
	
  void save_triangulation()
  {
    QString fileName = 
      QFileDialog::getSaveFileName( "triangulation.cgal", 
				  "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {
      // got a file name
      std::ofstream out(fileName);
      CGAL::set_ascii_mode(out);
      out << tr1 << std::endl;    
    }
  }

	

  void load_triangulation()
  {
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    tr1.clear();
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> tr1;
    something_changed();
  }

private:
  void show_conflicts(Point p)
  {
    if(tr1.dimension()<2) return;
    std::list<Face_handle> conflict_faces;
    std::list<Edge>  hole_bd;
    tr1.get_conflicts_and_boundary(p, 
    std::back_inserter(conflict_faces),
    std::back_inserter(hole_bd));
    std::list<Face_handle>::iterator fit = conflict_faces.begin();
    std::list<Edge>::iterator eit = hole_bd.begin();
    *widget << CGAL::WHITE ;
    for( ; fit != conflict_faces.end(); fit++)  {
      if(! tr1.is_infinite( *fit))
	*widget << tr1.triangle( *fit );
    }
    *widget << CGAL::YELLOW;
    for( ; eit != hole_bd.end(); eit++)  {
      if(! tr1.is_infinite( *eit ))
	*widget << tr1.segment( *eit );
    }		
  }
  inline  void something_changed(){current_state++;};


  CGAL::Qt_widget	                  *widget;		
  CGAL::Tools_toolbar	              *newtoolbar;
  CGAL::Layers_toolbar	            *vtoolbar;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  bool			  got_point;	
	  //if a CGAL::Point is received should be true
  int			  old_state;
};//endclass

#include "triangulation_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
  Window W(600,600); // physical widgetdow size
  app.setMainWidget(&W);
  W.setCaption(my_title_string);
  W.setMouseTracking(TRUE);
  W.show();
  // because Qt send resizeEvent only on show.
  W.set_window(-1, 1, -1, 1);
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
