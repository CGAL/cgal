// ============================================================================
//
// Copyright (c) 1997-2003 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
// ----------------------------------------------------------------------
//
// file          : polygon_2.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl; return 0;
}
#else

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include "polygon_2_toolbar.h"
#include <CGAL/IO/pixmaps/demoicon.xpm>
//#include <CGAL/IO/Qt_widget_get_polygon.h>


#include <qapplication.h>
#include <qmainwindow.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtextbrowser.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtimer.h>

#include <list>
#include <fstream>
#include <stack>
#include <set>
#include <string>


//global flags and variables
bool         is_point_visible = false;
int          current_state;
Cgal_Polygon polygon;
Point        point;

const QString my_title_string("Polygon_2 Demo with"
			      " CGAL Qt_widget");

class Qt_layer_show_polygon : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_polygon(){};
  void draw()
  {
    widget->lock();
      *widget << CGAL::LineWidth(3);
      *widget << CGAL::BLUE;
      *widget << polygon;
      *widget << CGAL::LineWidth(1);
      *widget << CGAL::WHITE;
      *widget << polygon;
      *widget << CGAL::GREEN;
      if(is_point_visible)
        *widget << point;
    widget->unlock();
  };
};//end class 


class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h){
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Polygon", this, SLOT(load_polygon()), CTRL+Key_L);
    file->insertItem("&Save Polygon", this, SLOT(save_polygon()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Edit", draw );
    draw->insertItem("Generate Polygon", this, SLOT(gen_poly()),
		     CTRL+Key_G);
    draw->insertItem("Show Polygon Information", this,
		     SLOT(show_info()), CTRL+Key_I);
    draw->insertItem("Show Point Information", this,
		     SLOT(show_pinfo()), CTRL+Key_F);

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the polygon_toolbar
    pt = new Polygon_toolbar(widget, this);

    *widget << CGAL::LineWidth(2) << CGAL::PointSize(6) 
	    <<CGAL::PointStyle(CGAL::DISC)
	    << CGAL::BackgroundColor (CGAL::BLACK);
  
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);

    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

    //layers
    widget->attach(&testlayer);

    qte = new QTextBrowser(NULL, "INFO");
    qte->setCaption("Information Window");
    //    qte->setTextFormat(PlainText);
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void new_instance()
  {
    widget->lock();
    stoolbar->clear_history();
    polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
		// set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:
  void gen_poly(){
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    // set the Visible Area to the Interval
    polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
        CGAL::random_polygon_2(100,
    			   std::back_inserter(polygon),
    			   Point_generator(1));
    show_info();
    something_changed();
  }

  void get_new_object(CGAL::Object obj)
  {
    Cgal_Polygon poly;
    Point p;
    if(CGAL::assign(poly, obj)) {
      polygon = poly;
      something_changed();
      show_info();
    } else if(CGAL::assign(p, obj)) {
      point = p;
      is_point_visible = true;
      show_pinfo();
      something_changed();
    }
  };


void show_pinfo(){
  qte->resize(300, 100);
  qte->show();
  qte->setText("Information on point");

  if(!polygon.is_simple()){
    qte->append("  The polygon is not simple!");
    return;
  }
  CGAL::Bounded_side bside   = polygon.bounded_side(point);
  switch (bside) {
    case CGAL::ON_BOUNDED_SIDE:
      qte->append("  The point is inside the polygon"); break;
    case CGAL::ON_BOUNDARY:
      qte->append("  The point is on the boundary of the polygon"); break;
    case CGAL::ON_UNBOUNDED_SIDE:
      qte->append("  The point is outside the polygon"); break;
  }
}

//--------------------------------------------------------------------------//
//                   PrintPolygonInfo
//--------------------------------------------------------------------------//
// prints some information about the polygon P to cerr

void show_info()
{
  QString s("%1");
  qte->resize(300, 300);
  qte->show();
  qte->setText("Polygon Information");

  if(polygon.is_empty())
    {qte->resize(300, 100); qte->append("P is empty!"); return;}
  s.setNum(polygon.size(), 10);
  qte->append("P.size() = " + s);
  s.setNum(polygon.is_empty(), 10);
  qte->append("P.is_empty() = " + s);
  s.setNum(polygon.is_simple(), 10);
  qte->append("P.is_simple() = " + s);
  s.setNum(polygon.is_convex(), 10);
  qte->append("P.is_convex() = " + s);

  if(polygon.is_simple()){
    CGAL::Orientation o = polygon.orientation();
    switch (o) {
      case CGAL::CLOCKWISE:
          qte->append("P.orientation() = CLOCKWISE");
	  break;
      case CGAL::COUNTERCLOCKWISE:
          qte->append("P.orientation() = COUNTERCLOCKWISE");
          break;
      case CGAL::COLLINEAR:
          qte->append("P.orientation() = COLLINEAR");
          break;
    }
  }
  QString s1, s2, s3, s4;
  s1.setNum(polygon.bbox().xmin());
  s2.setNum(polygon.bbox().ymin());
  s3.setNum(polygon.bbox().xmax());
  s4.setNum(polygon.bbox().ymax());
  s = "xmin:" + s1 + ", ymin:" + s2 + ", xmax:" + s3 + ", ymax:" + s4;
  qte->append("P.bbox(): " + s);
  s.setNum(float(polygon.area()));
  qte->append("P.area() = " + s);

  s1.setNum(float((*(polygon.left_vertex())).x()));
  s2.setNum(float((*(polygon.left_vertex())).y()));
  s = s1 + ", " + s2;
  qte->append("P.left_vertex() = " + s);
  s1.setNum(float((*(polygon.right_vertex())).x()));
  s2.setNum(float((*(polygon.right_vertex())).y()));
  s = s1 + ", " + s2;
  qte->append("P.right_vertex() = " + s);
  s1.setNum(float((*(polygon.top_vertex())).x()));
  s2.setNum(float((*(polygon.top_vertex())).y()));
  s = s1 + ", " + s2;
  qte->append("P.top_vertex() = " + s);
  s1.setNum(float((*(polygon.bottom_vertex())).x()));
  s2.setNum(float((*(polygon.bottom_vertex())).y()));
  s = s1 + ", " + s2;
  qte->append("P.bottom_vertex() = " + s);
}


  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Polygon\n"
  		"Copyright CGAL @2003");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
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

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void save_polygon()
  {
    QString fileName = QFileDialog::getSaveFileName( 
		"polygon.cgal", "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {                 // got a file name
      std::ofstream out(fileName);
      //out << std::setprecision(15);
      CGAL::set_ascii_mode(out);
      out << polygon << std::endl;
    }
  }

  void load_polygon()
  {
    QString s( QFileDialog::getOpenFileName(
		QString::null, "CGAL files (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    in >> polygon;
    something_changed();
  }

private:
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar
                         *stoolbar;
  Polygon_toolbar        *pt;
  int                    old_state;
  Qt_layer_show_polygon  testlayer;
  QTextBrowser           *qte;
};

#include "polygon_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon); 
  widget.show();
  current_state = -1;
  return app.exec();
}


#endif
