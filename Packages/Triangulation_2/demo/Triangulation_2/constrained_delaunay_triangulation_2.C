// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for 
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : constrained_delaunay_triangulation_2.C
// package       : Qt_widget
// author(s)     : Radu Ursu 
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


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
#include <list>
#include <set>
#include <string>

#include "constrained_cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "constrained_delaunay_triangulation_2_toolbar.h"
#include "constrained_delaunay_triangulation_2_toolbar_layers.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qprogressbar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>

const QString my_title_string("Constrained Delaunay Triangulation Demo with"
			      " CGAL Qt_widget");

CDT   ct;
int   current_state;
Coord_type        xmin, ymin, xmax, ymax;

class Window : public QMainWindow
{
  Q_OBJECT
public:
  Window(int w, int h)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
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
    file->insertItem("&Load Constraints", this, 
		      SLOT(load_constraints()), CTRL+Key_C);
    file->insertItem("&Save Constraints", this,
		      SLOT(save_constraints()), CTRL+Key_T);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()));
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
    newtoolbar = new Tools_toolbar(widget, this, &ct);	
    //the new scenes toolbar
    vtoolbar = new Layers_toolbar(widget, this, &ct);
  
    *widget << CGAL::BackgroundColor (CGAL::BLACK);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);

    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;
  };

  void  init_coordinates(){
    xmin = -1; xmax = 1;
    ymin = -1; ymax = 1;
  }


private slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    stoolbar->clear_history();
    ct.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
    something_changed();
  }
	
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    Segment s;
    Cgal_Polygon poly;
    if (CGAL::assign(s,obj))
    {
      Vertex_handle vs = ct.insert(s.source());
      Vertex_handle vt = ct.insert(s.target());
      ct.insert_constraint(vs, vt);
      something_changed();
    } else if(CGAL::assign(p,obj)) {
      ct.insert(p);
      something_changed();
    } else if (CGAL::assign(poly, obj)) {
      typedef Cgal_Polygon::Vertex_const_iterator VI;
      VI i=poly.vertices_begin();
      Point lp(i->x(), i->y()), lp1(lp);
      if(i!=poly.vertices_end()) {
        i++;
        for(;i!=poly.vertices_end();i++){
          Point p(i->x(), i->y());
          Vertex_handle vs = ct.insert(p);
          Vertex_handle vt = ct.insert(lp);
          ct.insert_constraint(vs, vt);
          lp = p;
        }
      }
      Vertex_handle vs = ct.insert(lp);
      Vertex_handle vt = ct.insert(lp1);
      ct.insert_constraint(vs, vt);
      something_changed();
    }
  }

  void howto(){
    QString home;
    home = "help/cindex.html";
    CGAL::Qt_help_window *help = new 
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Constrained Delaunay Triangulation,\n"
  		"Copyright CGAL @2002");
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    Window *ed = new Window(500, 500);
    ed->setCaption("Layer");    
    if(ct.number_of_vertices() > 1){
    Vertex_iterator it = ct.vertices_begin();
    xmin = xmax = (*it).point().x();
    ymin = ymax = (*it).point().y();
    while(it != ct.vertices_end()) {
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
    ed->stoolbar->clear_history();
    ed->widget->set_window(xmin, xmax, ymin, ymax);
    ed->show();
    }
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
    ct.clear();
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<200; count++)
      ct.insert(*g++);
    Vertex_iterator it = ct.vertices_begin();
    xmin = xmax = (*it).point().x();
    ymin = ymax = (*it).point().y();
    while(it != ct.vertices_end()) {
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
	
  void save_triangulation()
  {
    QString fileName = 
      QFileDialog::getSaveFileName( "triangulation.cgal", 
				  "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {
      // got a file name
      std::ofstream out(fileName);
      CGAL::set_ascii_mode(out);
      out << ct << std::endl;    
    }
  }

  void load_triangulation()
  {
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    ct.clear();
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> ct;
    Vertex_iterator it = ct.vertices_begin();
    xmin = xmax = (*it).point().x();
    ymin = ymax = (*it).point().y();
    while(it != ct.vertices_end()) {
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

  void load_constraints()
  {
    QString s( QFileDialog::getOpenFileName( 
	     QString::null,
	     "Constrained edges (*.edg);;Shewchuck Triangle .poly files (*.poly);;Constraints File (*.cst);;All files (*)", 
	     this ) );
    if ( s.isEmpty() ){
      return;
    }
    if(s.right(4) == ".edg"){
      std::ifstream in(s);
      CGAL::set_ascii_mode(in);

      QProgressBar *progress = new QProgressBar(NULL, "MyProgress");
      progress->setCaption("Loading constraints");
      int nedges = 0;
      ct.clear();
      in>>nedges;
      progress->setTotalSteps(nedges);
      progress->show();
      for(int n = 0; n<nedges; n++) {
        Point p1, p2;
        in >> p1 >> p2;
	      if(n==0){
	        xmin = xmax = p1.x();
	        ymin = ymax = p1.y();
	      }
	      if(xmin > p1.x())
	        xmin = p1.x();
	      if(xmax < p1.x())
	        xmax = p1.x();
	      if(ymin > p1.y())
	        ymin = p1.y();
	      if(ymax < p1.y())
	        ymax = p1.y();
	      if(xmin > p2.x())
	        xmin = p2.x();
	      if(xmax < p2.x())
	        xmax = p2.x();
	      if(ymin > p2.y())
	        ymin = p2.y();
	      if(ymax < p2.y())
	        ymax = p2.y();
        Vertex_handle vs = ct.insert(p1);
        Vertex_handle vt = ct.insert(p2);
        ct.insert_constraint(vs, vt);
        progress->setProgress(n);
      }
      progress->setProgress(nedges);
      progress->hide();
    } else if(s.right(5) == ".poly"){

    } else if(s.right(4) == ".cst"){
      std::ifstream in(s);
      CGAL::set_ascii_mode(in);
      std::istream_iterator<Point> it(in), done;
      bool first(true);
      CGAL::Bbox_2 b;
 
      Vertex_handle p_vh;
      Point p_q;

      while(it != done){
        Point p(*it);
        if(first){
          b = p.bbox();
        } else {
          b = b + p.bbox();
        }
        ++it;
        Point q(*it);
        b = b + q.bbox();
        ++it;
        if( (! first) && (p_q == p)){
          Vertex_handle vh = ct.insert(q);
          ct.insert_constraint(p_vh, vh);
          p_vh = vh;
        } else {
          Vertex_handle vhp = ct.insert(p);
          p_vh = ct.insert(q, vhp->face());
        }
        p_q = q;
        first = false;
      }
    }
    stoolbar->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
    something_changed();
  }

  void save_constraints()
  {
  }

private:
  inline  void something_changed(){current_state++;};


  CGAL::Qt_widget                   *widget;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  Tools_toolbar                     *newtoolbar;
  Layers_toolbar                    *vtoolbar;
  int                               old_state;
};//endclass

#include "constrained_delaunay_triangulation_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  Window W(500,500); // physical widgetdow size
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
