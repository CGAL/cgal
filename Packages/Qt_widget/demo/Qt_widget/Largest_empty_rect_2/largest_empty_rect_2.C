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
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : demo/Qt_widget/Max_k-gon/max_k-gon.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : CGAL-2.4
// release_date  : 2002, May 16
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

#include "cgal_types.h"

#include <fstream>
#include <stack>
#include <set>
#include <string>
#include <list>


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include "Qt_widget_toolbar.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_layer.h>

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



const QString my_title_string("Largest_empty_rectangle_2 Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int                 current_state;
std::list<Point>	  list_of_points;
Iso_rectangle_2     iso_rect;


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_ch(){};


  void draw()
  {
    widget->lock();
      //The ISO_RECTANGLE
      *widget << CGAL::RED << CGAL::LineWidth(3);
      *widget << CGAL::FillColor(CGAL::BLACK);
      *widget << iso_rect;
      *widget << CGAL::RED << CGAL::LineWidth(2);
      
      //The Largest_empty_rectangle
      Largest_empty_rect empty_rectangle(iso_rect);

      //VERTICES
      *widget << CGAL::PointSize(7) << CGAL::PointStyle(CGAL::CROSS);
      *widget << CGAL::GREEN;
      std::list<Point>::iterator itp = list_of_points.begin();
      while(itp!=list_of_points.end()){
        empty_rectangle.insert(*itp);
        *widget << (*itp++);
      }
      if(list_of_points.size() > 0) {
        //Show the Largest_emty_rectangle
        Iso_rectangle_2 b = empty_rectangle.get_largest_empty_iso_rectangle();
        *widget << CGAL::BLUE << CGAL::LineWidth(1) 
              << CGAL::FillColor(CGAL::BLUE);
        *widget << b;
      }

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
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate points", this,
				SLOT(gen_points()), CTRL+Key_G );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the new tools toolbar
    //setUsesBigPixmaps(TRUE);
    newtoolbar = new CGAL::Tools_toolbar(widget, this, &list_of_points);	
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
  
    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::GRAY);
  
    resize(w,h);
    widget->show();

    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
    this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

    //layers
    widget->attach(&testlayer);

    //initialize the iso_rectangle
    iso_rect = Iso_rectangle_2(Point(-0.8, -0.8), Point(0.8, 0.8));
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  void new_instance()
  {
    widget->lock();
    list_of_points.clear();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
			// set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point p;
    if(CGAL::assign(p,obj)) {
      list_of_points.push_back(p);
      something_changed();
    }
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Maximum inscribed k-gon,\n"
  		"Copyright CGAL @2002");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->show();
    ed->set_window(-1.1, 1.1, -1.1, 1.1);
    something_changed();
  }

  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }	

  void gen_points()
  {
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
		// set the Visible Area to the Interval

    // send resizeEvent only on show.
    CGAL::Random_points_in_disc_2<Point> g(0.5);
    for(int count=0; count<100; count++) {
      list_of_points.push_back(*g++);
    }
    something_changed();
  }
	
	

private:
  CGAL::Qt_widget       *widget;		
  CGAL::Tools_toolbar   *newtoolbar;
  CGAL::Qt_widget_standard_toolbar
                        *stoolbar;
  int                   old_state;
  Qt_layer_show_ch      testlayer;
};

#include "largest_empty_rect_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  widget.show();
  // because Qt send resizeEvent only on show.
  widget.set_window(-1, 1, -1, 1);
  current_state = -1;
  return app.exec();
}

#endif // CGAL_USE_QT
