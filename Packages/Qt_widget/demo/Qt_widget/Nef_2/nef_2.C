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
// file          : demo/Qt_widget/Nef_2/nef_2.C
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

#include <fstream>
#include <stack>
#include <set>
#include <list>

#include "cgal_types.h"
#include "nef_2.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Nef_2.h>
#include "Qt_widget_toolbar.h"

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


const QString my_title_string("Nef_2 Demo with"
			      " CGAL Qt_widget");

//global flags and variables
int current_state;
Nef_polyhedron Nef_visible(Nef_polyhedron::EMPTY);
Nef_polyhedron Nef_visible2(Nef_polyhedron::EMPTY);


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_ch(){};


  void draw()
  {
    widget->setRasterOp(XorROP);
	  *widget << CGAL::FillColor(CGAL::GRAY) << CGAL::WHITE;
    *widget << Nef_visible;
    *widget << CGAL::FillColor(CGAL::RED) << CGAL::GREEN;
    *widget << Nef_visible2;
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
    file->insertItem("Load Nef_2", this, SLOT(load_nef()), CTRL+Key_L);
    file->insertItem("Save Nef_2", this, SLOT(save_nef()), CTRL+Key_S);
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
    newtoolbar = new CGAL::Tools_toolbar(widget, this);	
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);
    this->addToolBar(newtoolbar->toolbar(), Top, FALSE);
  
    Line l1(Point(0, 0), Point(0, 2));
    Nef_polyhedron N1(l1, Nef_polyhedron::INCLUDED);
    Nef_visible = N1;

    Line l2(Point(0, 0), Point(2, 0));
    Nef_polyhedron N2(l2, Nef_polyhedron::INCLUDED);
    Nef_visible2 = N2;

    *widget << CGAL::BackgroundColor (CGAL::BLACK);
  
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
    widget->set_window(-1.1, 1.1, -1.1, 1.1); 
			// set the Visible Area to the Interval
	Nef_polyhedron N_temp(Nef_polyhedron::EMPTY);
	Nef_visible = N_temp;
    widget->unlock();
    something_changed();
  }


  void load_nef()
	{
		QString s( QFileDialog::getOpenFileName( QString::null,
		    "CGAL files (*.cgal)", this ) );
		if ( s.isEmpty() )
			return;
		std::ifstream in(s);
		CGAL::set_ascii_mode(in);
		
		Nef_polyhedron N_temp(Nef_polyhedron::EMPTY);
		Nef_visible = N_temp;
		in >> Nef_visible;
		something_changed();
	}

  void save_nef()
	{
		QString fileName = 
		QFileDialog::getSaveFileName( "nef_2.cgal", 
				  "Cgal files (*.cgal)", this );
		if ( !fileName.isNull() ) {
			// got a file name
			std::ofstream out(fileName);
			CGAL::set_ascii_mode(out);
			out << Nef_visible << std::endl;    
		}
	}//end save_nef()

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Cartesian_point_2   p;
    Cartesian_polygon_2 poly;
    if(CGAL::assign(p, obj)) {
      CGAL::Quotient<RT> wsxq = double_to_quotient<RT>(p.x());
      CGAL::Quotient<RT> wsyq = double_to_quotient<RT>(p.y());
      RT wsx = wsxq.numerator() * wsyq.denominator(); 
      RT wsy = wsyq.numerator() * wsxq.denominator(); 
      RT wsh  = wsxq.denominator() * wsyq.denominator(); 
      Point p1(wsx, wsy, wsh);
      Point pt[1] = {p1};
      Nef_polyhedron Nt(pt, pt+1);
      Nef_visible = Nef_visible + Nt;
      something_changed();
    } else if(CGAL::assign(poly, obj)){
      Vertex_iterator it = poly.vertices_begin();
      std::list<Point> l_of_p;
      while(it != poly.vertices_end()){
        double xp = (*it).x();
        double yp = (*it).y();
        CGAL::Quotient<RT> wsxq = double_to_quotient<RT>(xp);
        CGAL::Quotient<RT> wsyq = double_to_quotient<RT>(yp);
        RT wsx = wsxq.numerator() * wsyq.denominator(); 
        RT wsy = wsyq.numerator() * wsxq.denominator(); 
        RT wsh  = wsxq.denominator() * wsyq.denominator(); 
        Point p1(wsx, wsy, wsh);
        l_of_p.push_back(p1);
        it++;
      }
      Nef_polyhedron Nt(l_of_p.begin(), l_of_p.end(), Nef_polyhedron::INCLUDED);
      Nef_visible = Nef_visible + Nt;
    }
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Nef_2,\n"
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

private:
  CGAL::Qt_widget			*widget;		
  CGAL::Tools_toolbar		*newtoolbar;
  CGAL::Qt_widget_standard_toolbar  *stoolbar;
  int						old_state;  	
  Qt_layer_show_ch			testlayer;

};

#include "nef_2.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
    app.setStyle( new QPlatinumStyle );
    QPalette p( QColor( 250, 215, 100 ) );
    app.setPalette( p, TRUE );
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
