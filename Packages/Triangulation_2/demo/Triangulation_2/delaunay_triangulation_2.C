// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : delaunay_triangulation_2.C
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

//Application headers
#include "cgal_types.h"
#include "delaunay_triangulation_2_toolbar.h"
#include "delaunay_triangulation_2_toolbar_layers.h"

//Qt_widget headers
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

//STL headers
#include <fstream>
#include <stack>
#include <set>
#include <string>

//Qt headers
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

const QString my_title_string("Delaunay Triangulation Demo with"
			      " CGAL Qt_widget");

Delaunay	tr1;
int		current_state;
Coord_type      xmin, ymin, xmax, ymax;

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
		      SLOT(save_triangulation()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp,
		      SLOT( closeAllWindows() ), CTRL+Key_Q );


    // edit menu
    QPopupMenu * edit = new QPopupMenu( this );
    menuBar()->insertItem( "&Edit", edit );
    edit->insertItem("&Generate_triangulation", this,
		      SLOT(generate_triangulation()), CTRL+Key_G );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this, &tr1);	
    //the new scenes toolbar
    vtoolbar = new Layers_toolbar(widget, this, &tr1);
  
    *widget << CGAL::BackgroundColor (CGAL::BLACK);
    *widget << CGAL::LineWidth(2);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);

    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    connect(newtoolbar, SIGNAL(changed()), 
	    this, SLOT(something_changed()));

    //application flag stuff
    got_point = FALSE;
    old_state = 0;
    triangulation_changed = true;
  };

  void  init_coordinates(){
    xmin = -1; xmax = 1;
    ymin = -1; ymax = 1;
  }

private slots:
  void new_instance(){
    widget->lock();
    widget->clear();
    stoolbar->clear_history();
    tr1.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
    triangulation_changed = true;
    something_changed();
  }
	
  void get_new_object(CGAL::Object obj){
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
      triangulation_changed = true;
    } 
  }

  void insert_after_show_conflicts(QMouseEvent*){
    if(got_point)
    {
      got_point = FALSE;
      something_changed();
    }
  }

  void howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new 
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void about(){
    QMessageBox::about( this, my_title_string,
		"This is a demo for Delaunay Triangulation 2,\n"
  		"Copyright CGAL @2001");
  }

  void aboutQt(){
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    Window *ed = new Window(500, 500);
    ed->setCaption("Layer");    
    if(tr1.number_of_vertices() > 1){
      Vertex_iterator it = tr1.vertices_begin();
      xmin = xmax = (*it).point().x();
      ymin = ymax = (*it).point().y();
      while(it != tr1.vertices_end()) {
        if(xmin > (*it).point().x())
          xmin = (*it).point().x();
        if(xmax < (*it).point().x())
          xmax = (*it).point().x();
        if(ymin > (*it).point().y())
          ymin = (*it).point().y();
        if(ymax < (*it).point().y())
          ymax = (*it).point().y();
        it++;
      }
    }
    ed->stoolbar->clear_history();
    ed->widget->set_window(xmin, xmax, ymin, ymax);
    ed->show();
    something_changed();
  }

  void timerDone(){
    if(triangulation_changed){
      if(tr1.number_of_vertices() > 2)
        newtoolbar->set_line_enabled(true);
      else
        newtoolbar->set_line_enabled(false);
      if(tr1.number_of_vertices() > 2)
        newtoolbar->set_move_enabled(true);
      else
	newtoolbar->set_move_enabled(false);
      triangulation_changed = false;
    }
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
      triangulation_changed = true;
    }
  }	

  void generate_triangulation(){
    tr1.clear();
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++)
      tr1.insert(*g++);
    Vertex_iterator it = tr1.vertices_begin();
    xmin = xmax = (*it).point().x();
    ymin = ymax = (*it).point().y();
    while(it != tr1.vertices_end()) {
      if(xmin > (*it).point().x())
        xmin = (*it).point().x();
      if(xmax < (*it).point().x())
        xmax = (*it).point().x();
      if(ymin > (*it).point().y())
        ymin = (*it).point().y();
      if(ymax < (*it).point().y())
        ymax = (*it).point().y();
      it++;
    }
    stoolbar->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
    triangulation_changed = true;
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

    Vertex_iterator it = tr1.vertices_begin();
    xmin = xmax = (*it).point().x();
    ymin = ymax = (*it).point().y();
    while(it != tr1.vertices_end()) {
      if(xmin > (*it).point().x())
        xmin = (*it).point().x();
      if(xmax < (*it).point().x())
        xmax = (*it).point().x();
      if(ymin > (*it).point().y())
        ymin = (*it).point().y();
      if(ymax < (*it).point().y())
        ymax = (*it).point().y();
      it++;
    }
    stoolbar->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
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

public slots:
  inline  void something_changed(){current_state++;};

private:
  CGAL::Qt_widget                   *widget;		
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  Tools_toolbar                     *newtoolbar;
  Layers_toolbar                    *vtoolbar;
  bool                              got_point;	
  //if a CGAL::Point is received should be true
  bool                              triangulation_changed;
  //true only when triangulation has changed
  int                               old_state;
};//endclass

#include "delaunay_triangulation_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  Window W(600,600); // physical widgetdow size
  app.setMainWidget(&W);
  W.setCaption(my_title_string);
  W.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
  W.show();  
  W.init_coordinates();  
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
